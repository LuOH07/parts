#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
质粒 .dna → iGEM 砖号 txt
=========================

输入质粒文件（SnapGene / SnapGene Viewer 保存的 .dna），输出 .txt：

    一行一段，按质粒自身的坐标顺序排列（从 1 号碱基走到末尾），每行是

        名字（正向）：砖号          这一段查到了砖，而且是正链
        名字（反向）：DNA 序列      这一段查不到砖号，直接写序列

    所以整条质粒的每一段都会出现在输出里，不会漏。
    没被标注覆盖的区域，名字写成「无标注区(起..止)」。
    标「反向」的段，写出来的是它自己的 5'→3'（也就是质粒正链的反向互补）。

可以边查边加：一个质粒正在查的时候，直接再选下一个，
它会自己再开一个窗口同时查（最多 3 个，多的排队）。

· 读 .dna 的逻辑来自 snapgene-to-txt.html。
· 查砖号的逻辑来自 igem_brick_finder.py（沿用 app/login-state.json 的登录态）。
"""

from __future__ import annotations

import argparse
import gzip
import html as _html
import os
import queue
import re
import shutil
import struct
import sys
import threading
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

# 查砖号（登录态、BLAST 流程）沿用原来的程序
import igem_brick_finder as igem  # noqa: E402

APP_VERSION = "2.0"

MIN_QUERY_LEN = 7          # BLAST 接口要求至少 7 bp
MAX_QUERY_LEN = 100000     # BLAST 接口上限
DEFAULT_JOBS = 3           # 默认同时查几个质粒

Event = Callable[..., None]


# ==========================================================================
# 1. SnapGene .dna 解析（从 snapgene-to-txt.html 移植）
# ==========================================================================

SEG_COOKIE = 0x09
SEG_DNA = 0x00
SEG_PRIMERS = 0x05
SEG_NOTES = 0x06
SEG_FEATURES = 0x0A

_RANGE_RE = re.compile(r"^(-?\d+)-(-?\d+)$")
_GC_RE = re.compile(r"\s+")


@dataclass
class Feature:
    """一条标注（Feature 段里的 <Feature>，或 Primers 段里的一个结合位点）。"""

    type: str
    name: str
    strand: int                      # 1 正向 / -1 反向
    segs: List[Tuple[int, int]]      # 1-based，含首含尾
    swapped: bool = False
    source: str = "feature"          # feature / primer
    order: int = 0

    @property
    def label(self) -> str:
        return self.name.strip() if self.name and self.name.strip() else self.type


@dataclass
class Parsed:
    seq: str
    topology: str                    # circular / linear
    features: List[Feature]
    feature_count: int = 0
    primer_count: int = 0


def read_segments(data: bytes) -> List[Tuple[int, bytes]]:
    """把文件读成一串 [type:1][len:4 大端][payload]。"""
    segs: List[Tuple[int, bytes]] = []
    i = 0
    n = len(data)
    while i + 5 <= n:
        seg_type = data[i]
        length = struct.unpack(">I", data[i + 1:i + 5])[0]
        i += 5
        if i + length > n:           # 截断的坏段：停止
            break
        segs.append((seg_type, data[i:i + length]))
        i += length
    return segs


def _parse_xmlish(payload: bytes) -> List[str]:
    """把可能被二进制/压缩块包裹的 XML 抠出来。"""
    texts: List[str] = []
    n = len(payload)
    i = 0
    while i < n:
        b = payload[i]
        if b < 0x20 and b not in (0x09, 0x0A, 0x0D):
            i += 1
            continue
        j = i
        while j < n:
            b = payload[j]
            if b < 0x20 and b not in (0x09, 0x0A, 0x0D):
                break
            j += 1
        if j - i >= 20:
            chunk = payload[i:j]
            try:
                texts.append(chunk.decode("utf-8"))
            except UnicodeDecodeError:
                texts.append(chunk.decode("latin-1"))
        i = j

    out: List[str] = []
    for text in texts:
        k = text.find("<Features")
        if k == -1:
            k = text.find("<Primers")
        if k == -1:
            k = text.find("S{")
        if k == -1:                   # 不含 XML 元素的块（如 Descriptions）
            continue
        out.append(text[k:])
    return out


def _is_seq_byte(c: int) -> bool:
    return (c in (0x0A, 0x0D, 0x20, 0x2D, 0x2E)
            or 0x41 <= c <= 0x5A or 0x61 <= c <= 0x7A)


def _clean_seq(text: str) -> str:
    return re.sub(r"[^A-Za-z]", "", text).upper()


def _extract_seq_text(payload: bytes, seq_len_hint: int) -> str:
    """从序列段 payload 里取出序列本体（取最长的连续序列字符块）。"""
    n = len(payload)
    best = ""
    i = 0
    while i < n:
        if not _is_seq_byte(payload[i]):
            i += 1
            continue
        j = i
        while j < n and _is_seq_byte(payload[j]):
            j += 1
        if j - i > len(best):
            best = payload[i:j].decode("latin-1")
        i = j

    seq = _clean_seq(best)
    hint_ok = seq_len_hint >= 10 and (
        seq_len_hint + 1 == len(payload)
        or seq_len_hint * 2 + 1 == len(payload)
        or seq_len_hint + 1 == max(0, len(payload) - 1)
        or seq_len_hint * 2 + 1 == len(payload) - 1
    )
    if hint_ok and len(seq) >= seq_len_hint:
        # 环状序列在文件里存了两遍（第二遍开头补 N），按 cookie 长度裁掉多余部分
        seq = seq[:seq_len_hint]
    return seq


def _xunesc(s: str) -> str:
    return _html.unescape(s)


def _attr(tag: str, name: str) -> Optional[str]:
    m = re.search(name + r'\s*=\s*"([^"]*)"', tag)
    return _xunesc(m.group(1)) if m else None


def _parse_range(text: Optional[str]) -> Optional[Tuple[int, int]]:
    if not text:
        return None
    m = _RANGE_RE.match(_GC_RE.sub("", text))
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def _parse_feature(src: str) -> Optional[Feature]:
    open_m = re.search(r"<Feature\b[^>]*>", src)
    if not open_m:
        return None
    open_tag = open_m.group(0)
    direction = _attr(open_tag, "directionality") or "1"
    try:
        strand = -1 if int(direction) == 2 else 1
    except ValueError:
        strand = 1

    segs: List[Tuple[int, int]] = []
    sub_labels: List[str] = []
    for m in re.finditer(r"<Segment\b([^>]*?)/?>", src):
        tag = m.group(0)
        if (_attr(tag, "type") or "standard") == "gap":
            continue
        rng = _parse_range(_attr(tag, "range"))
        if rng is None:
            continue
        sub = _attr(tag, "name")
        if sub:
            sub_labels.append(sub)
        segs.append(rng)
    if not segs:
        return None

    label = _attr(open_tag, "label")
    if not label:
        label = _attr(open_tag, "name")
    if not label:
        label = next((s for s in sub_labels if s and s.strip()), None)

    return Feature(
        type=_attr(open_tag, "type") or "misc_feature",
        name=_GC_RE.sub(" ", label or "").strip(),
        strand=strand,
        segs=segs,
        swapped=_attr(open_tag, "swappedSegmentNumbering") == "1",
    )


def _parse_primers(blocks: Sequence[str]) -> List[Feature]:
    feats: List[Feature] = []
    seen = set()
    for block in blocks:
        for pm in re.finditer(r"<Primer\b[^>]*>[\s\S]*?</Primer>", block):
            chunk = pm.group(0)
            open_tag = re.search(r"<Primer\b[^>]*>", chunk).group(0)
            name = _GC_RE.sub(" ", (_attr(open_tag, "name") or "")
                              .replace("<", "").replace(">", "")).strip()
            for sm in re.finditer(r"<BindingSite\b[^>]*?/?>", chunk):
                tag = sm.group(0)
                rng = _parse_range(_attr(tag, "location"))
                if rng is None:
                    continue
                # 引物坐标比特征坐标偏移 1
                s, e = rng[0] + 1, rng[1] + 1
                strand = -1 if (_attr(tag, "boundStrand") or "0") == "1" else 1
                key = (name, strand, s, e)
                if key in seen:
                    continue
                seen.add(key)
                feats.append(Feature(type="primer_bind", name=name, strand=strand,
                                     segs=[(s, e)], source="primer"))
    # 同一引物有多个位置时加 (2)(3)… 后缀
    counts: Dict[str, int] = {}
    for f in feats:
        counts[f.name] = counts.get(f.name, 0) + 1
    used: Dict[str, int] = {}
    for f in feats:
        if counts.get(f.name, 0) > 1:
            used[f.name] = used.get(f.name, 0) + 1
            f.name = f"{f.name} ({used[f.name]})"
    return feats


def _detect_topology(segs, data: bytes) -> str:
    for seg_type, payload in segs:
        if seg_type == SEG_DNA and payload:
            return "circular" if (payload[0] & 0x01) else "linear"
    head = data[:min(len(data), 40000)].decode("latin-1")
    return "circular" if re.search(r"\bcircular\b", head, re.IGNORECASE) else "linear"


def parse_snapgene(data: bytes, include_primers: bool = False) -> Parsed:
    segs = read_segments(data)
    if not segs:
        raise ValueError("文件为空，或者不是有效的 SnapGene .dna 文件（没解析到任何数据段）。")

    header_type, header_payload = segs[0]
    header_text = header_payload[:16].decode("latin-1")
    if not header_text.startswith("SnapGene") and header_type != SEG_COOKIE:
        raise ValueError(
            '这不是 SnapGene 的 .dna 文件（文件头不是 "SnapGene"）。\n'
            "如果你手上的是 .gb / .fasta 文本，请先用 SnapGene 另存为 .dna。"
        )

    seq_len_hint = 0
    if len(header_text) >= 10 and len(header_payload) >= 10:
        seq_len_hint = (header_payload[8] << 8) | header_payload[9]

    seq = ""
    for seg_type, payload in segs:
        if seg_type != SEG_DNA or len(payload) < 2:
            continue
        candidate = _extract_seq_text(payload, seq_len_hint)
        if len(candidate) > len(seq):
            seq = candidate
    if not seq:
        raise ValueError("文件里没有找到 DNA 序列段，无法继续。")

    raw_features: List[Feature] = []
    raw_primers: List[Feature] = []
    for seg_type, payload in segs:
        if seg_type not in (SEG_FEATURES, SEG_PRIMERS):
            continue
        for block in _parse_xmlish(payload):
            if seg_type == SEG_FEATURES:
                for fm in re.finditer(r"<Feature\b[^>]*>[\s\S]*?</Feature>", block):
                    f = _parse_feature(fm.group(0))
                    if f is not None:
                        raw_features.append(f)
            else:
                raw_primers.extend(_parse_primers([block]))

    for i, f in enumerate(raw_features):
        f.order = i
    features = list(raw_features)
    if include_primers:
        for i, f in enumerate(raw_primers):
            f.order = len(raw_features) + i
        features += raw_primers

    return Parsed(seq=seq, topology=_detect_topology(segs, data), features=features,
                  feature_count=len(raw_features), primer_count=len(raw_primers))


# ==========================================================================
# 2. 序列工具 + 把特征摊成"一行一段"的片段表
# ==========================================================================

_COMPLEMENT = str.maketrans("ATGCURYKMSWBDHVN", "TACGAYRMKSWVHDBN")


def revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT)[::-1]


def _ranges_len(ranges: Sequence[Tuple[int, int]]) -> int:
    return sum(max(0, b - a) for a, b in ranges)


def feature_ranges(f: Feature, length: int, circular: bool) -> List[Tuple[int, int]]:
    """把 1-based 的段转成 0-based 半开区间，按"阅读方向"排列。"""
    pieces: List[Tuple[Tuple[int, int], bool]] = []
    wrapped = False
    for raw_s, raw_e in f.segs:
        s, e = raw_s, raw_e
        if s < 1:
            s = 1
        # 跨起点：SnapGene 记作 start > end（320-30 表示 …320,321,1…30）
        if circular and raw_s > raw_e:
            # 但 "2135-2134" 这种 start = end+1 的写法是退化/空区间，
            # 不是跨起点：照跨起点算会把整条质粒都占掉，把别的标注全吃掉。
            e_eff = max(0, min(raw_e, length))
            wrapped_len = (length - s + 1) + e_eff
            if wrapped_len >= length:
                continue
            wrapped = True
            if s - 1 < length:
                pieces.append(((min(s - 1, length), length), True))
            if e >= 1:
                pieces.append(((0, min(e, length)), True))
            continue
        if e < s:
            e = s
        a = min(s - 1, max(0, length - 1))
        b = min(e, length)
        if a < b:
            pieces.append(((a, b), False))

    if not pieces:
        gs, ge = f.segs[0]
        a = max(0, min((gs or 1) - 1, length - 1))
        return [(a, min(length, max(a + 1, ge or a + 1)))]

    ranges = [p[0] for p in pieces]
    # 反向链的多段特征：文件里段通常按坐标升序排列，阅读顺序应从高坐标读起
    contiguous_asc = len(pieces) > 1 and all(
        i == 0 or pieces[i - 1][0][1] == pieces[i][0][0] for i in range(len(pieces)))
    if f.strand < 0 and len(f.segs) > 1 and not f.swapped and not wrapped and contiguous_asc:
        ranges = ranges[::-1]
    if len(ranges) > 1 and _ranges_len(ranges) > length:
        ranges = ranges[:-1]
    return ranges


def part_seq(seq: str, ranges: Sequence[Tuple[int, int]], strand: int) -> str:
    joined = "".join(seq[max(0, a):max(0, min(len(seq), b))] for a, b in ranges)
    return revcomp(joined) if strand < 0 else joined


def is_wrapped(f: Feature, circular: bool) -> bool:
    return bool(circular and len(f.segs) == 1 and f.segs[0][0] > f.segs[0][1])


@dataclass
class Piece:
    """输出里的一行。"""

    kind: str                        # "feature" / "gap"
    name: str
    ftype: str
    strand: int
    ranges: List[Tuple[int, int]]
    seq: str
    pos: int                         # 0-based，用来排序
    wrapped: bool = False
    brick: Optional[str] = None      # 查到的砖号

    @property
    def coords(self) -> str:
        if not self.ranges:
            return ""
        by_pos = sorted(self.ranges)
        if len(by_pos) == 1:
            return f"{by_pos[0][0] + 1}..{by_pos[0][1]}"
        if self.wrapped:
            return "跨起点 " + " + ".join(f"{a + 1}..{b}" for a, b in by_pos)
        return f"{by_pos[0][0] + 1}..{by_pos[-1][1]}"

    @property
    def strand_text(self) -> str:
        """这一行里的序列是正向还是反向。"""
        return "反向" if self.strand < 0 else "正向"

    @property
    def title(self) -> str:
        """输出行冒号前面的那部分：名字（正向/反向）。"""
        if self.kind == "feature":
            name = (self.name or "").strip()
            name = name or self.ftype or "未命名标注"
        else:
            name = f"无标注区({self.coords})"
        return f"{name}（{self.strand_text}）"


def _container_flags(groups: Sequence[Tuple[Feature, List[Tuple[int, int]]]],
                     length: int) -> set:
    """找出"只是给一整片区域起名字"的标注：整段包住了 ≥2 条别的标注。

    这种标注（比如把整个合成片段框起来的 putative synthesis）不是部件，
    要是让它参与切分，会把里面真正的部件全吃掉。
    """
    if len(groups) < 3 or length <= 0:
        return set()
    flags = set()
    for i, (_f, ranges) in enumerate(groups):
        own = _ranges_len(ranges)
        mark = bytearray(length)
        for a, b in ranges:
            mark[a:b] = b"\x01" * (b - a)
        # 前缀和，便于 O(1) 数区间里被盖住了几个碱基
        prefix = [0] * (length + 1)
        run = 0
        for x in range(length):
            run += mark[x]
            prefix[x + 1] = run
        inside = 0
        for j, (_g, granges) in enumerate(groups):
            if i == j:
                continue
            total = _ranges_len(granges)
            # 必须"比自己小"才算被包住；范围一样的重复标注不算
            if total and total < own and \
                    sum(prefix[b] - prefix[a] for a, b in granges) == total:
                inside += 1
                if inside >= 2:
                    flags.add(i)
                    break
    # 兜底：全都判成"区域"说明判断失效了，那就一条都不跳
    if len(flags) == len(groups):
        return set()
    return flags


def build_pieces(parsed: Parsed, *, rev_comp: bool = True,
                 log: Optional[Callable[[str], None]] = None) -> List[Piece]:
    """把特征摊成互不重叠的片段，加上没被覆盖的区域，按质粒坐标顺序排好。"""
    say = log or (lambda _m: None)
    seq = parsed.seq
    length = len(seq)
    circular = parsed.topology == "circular"
    covered = bytearray(length)

    candidates = []
    for f in parsed.features:
        # source 是「这条序列哪来的」，不是部件；它常常横跨或压住真部件
        if f.type == "source":
            say(f"跳过一条 source 标注（序列来源，不算部件）："
                f"{f.label} {','.join(f'{a}-{b}' for a, b in f.segs)}")
            continue
        ranges = feature_ranges(f, length, circular)
        if not ranges:
            continue
        candidates.append((f, ranges))

    # 整段包住别的标注的「区域标注」先剔掉，否则里面的部件全被它吞了
    containers = _container_flags(candidates, length)
    if containers:
        kept_candidates = []
        for idx, (f, ranges) in enumerate(candidates):
            if idx in containers:
                say(f"跳过一条「区域」标注（把别的标注整段包住了，不是部件）：{f.label}")
                continue
            kept_candidates.append((f, ranges))
        candidates = kept_candidates

    # 先按起点、再按长度（长的优先）、最后按原顺序
    prepared = []
    for f, ranges in candidates:
        prepared.append((min(a for a, _ in ranges), -_ranges_len(ranges), f.order, f, ranges))
    prepared.sort(key=lambda item: (item[0], item[1], item[2]))

    pieces: List[Piece] = []
    for _, _, _, f, ranges in prepared:
        # 裁掉已经被前面片段占掉的碱基：保证每条碱基只出现一次，且什么都不会丢。
        # 注意：这里必须保持 feature_ranges 给出的顺序（跨起点 / 反向链的段
        # 不是按坐标升序排的），按坐标排序会把序列读反。
        kept: List[Tuple[int, int]] = []
        for a, b in ranges:
            i = a
            while i < b:
                if covered[i]:
                    i += 1
                    continue
                j = i
                while j < b and not covered[j]:
                    j += 1
                kept.append((i, j))
                i = j
        if not kept:
            continue                     # 完全被前面的片段覆盖
        for a, b in kept:
            covered[a:b] = b"\x01" * (b - a)
        # 用 --no-revcomp 时输出的其实是正链，标记也要跟着改，别写反
        shown_strand = f.strand if rev_comp else 1
        text = part_seq(seq, kept, shown_strand)
        if not text:
            continue
        pieces.append(Piece(
            kind="feature", name=f.label, ftype=f.type, strand=shown_strand,
            ranges=kept, seq=text,
            pos=min(a for a, _ in kept),
            wrapped=is_wrapped(f, circular),
        ))

    # 没被任何特征覆盖的区域 → 直接输出序列
    i = 0
    while i < length:
        if covered[i]:
            i += 1
            continue
        j = i
        while j < length and not covered[j]:
            j += 1
        pieces.append(Piece(kind="gap", name="", ftype="", strand=1,
                            ranges=[(i, j)], seq=seq[i:j], pos=i))
        i = j

    pieces.sort(key=lambda p: p.pos)
    return pieces


def output_lines(pieces: Sequence[Piece]) -> List[str]:
    """一行一段：名字：砖号（查得到）或 名字：序列（查不到）。"""
    return [f"{p.title}：{p.brick or p.seq}" for p in pieces]


# ==========================================================================
# 3. 查砖号（沿用 igem_brick_finder 的登录态与 BLAST 流程）
# ==========================================================================

class LookupFailed(RuntimeError):
    pass


class SharedCache:
    """多个质粒之间共用已查过的结果（同一个砖常在不同质粒里重复出现）。"""

    MISS = object()

    def __init__(self) -> None:
        self._data: Dict[str, Optional[str]] = {}
        self._lock = threading.Lock()

    def get(self, key: str):
        with self._lock:
            return self._data.get(key, SharedCache.MISS)

    def put(self, key: str, value: Optional[str]) -> None:
        with self._lock:
            self._data[key] = value

    def __len__(self) -> int:
        with self._lock:
            return len(self._data)


def _pick_exact(hits: Sequence[igem.Hit], query_length: int) -> Optional[igem.Hit]:
    """在结果里找 100% 且等长的砖。

    站点把相似度排在前面的可能是 99.x% 的近亲（甚至只差一个 gap），
    所以要整表扫一遍，不能只看第一条。
    """
    best: Optional[igem.Hit] = None
    for hit in hits:
        if not igem.is_strict_exact(hit, query_length):
            continue
        if not hit.part_name:
            continue
        if best is None:
            best = hit
        elif (hit.bit_score, hit.identity) > (best.bit_score, best.identity):
            best = hit
    return best


class BrickFinder:
    """一个浏览器会话，连续查很多条序列。"""

    def __init__(self, *, profile_dir: Optional[Path] = None, headless: bool = True,
                 verbose: bool = False, timeout: float = 180.0,
                 progress: Optional[Callable[[str], None]] = None,
                 cache: Optional[SharedCache] = None,
                 allow_interactive_login: bool = True) -> None:
        self.profile_dir = Path(profile_dir) if profile_dir else igem.DEFAULT_PROFILE
        self.headless = headless
        self.verbose = verbose
        self.timeout = timeout
        self.progress = progress
        self.cache = cache
        self.allow_interactive_login = allow_interactive_login
        self.session: Optional[igem.Session] = None

    # -- 生命周期 ---------------------------------------------------------

    def _make_session(self, headless: bool) -> igem.Session:
        return igem.Session(
            headless=headless,
            profile_dir=self.profile_dir,
            timeout=self.timeout,
            verbose=self.verbose,
            progress=self.progress,
        )

    def set_progress(self, progress: Optional[Callable[[str], None]]) -> None:
        self.progress = progress
        if self.session is not None:
            self.session.progress = progress

    def open(self) -> None:
        """启动浏览器并确认登录。"""
        self.session = self._make_session(self.headless)
        self.session.start()
        try:
            self.session.open_search_page()
            self.session.ensure_login(interactive=False)
        except igem.LoginRequired:
            self.session.close()
            if not self.allow_interactive_login:
                raise
            # 登录态过期了：换一个能看见的窗口让用户登录一次
            self._progress("登录态已过期，打开浏览器窗口等你登录…")
            self.session = self._make_session(headless=False)
            self.session.start()
            self.session.open_search_page()
            self.session.ensure_login(interactive=True)

    def close(self) -> None:
        if self.session is not None:
            self.session.close()
        self.session = None

    def __enter__(self) -> "BrickFinder":
        try:
            self.open()
        except BaseException:
            # 打开失败也要把浏览器收掉，别留孤儿进程
            self.close()
            raise
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def _progress(self, msg: str) -> None:
        if self.progress is not None:
            try:
                self.progress(msg)
            except Exception:  # noqa: BLE001
                pass
        if self.verbose:
            print(f"[brick] {msg}", file=sys.stderr, flush=True)

    # -- 单条查询 ---------------------------------------------------------

    def query(self, sequence: str) -> Optional[str]:
        """查一条序列，返回砖号；查不到返回 None。"""
        if self.session is None:
            raise LookupFailed("浏览器会话没有启动。")
        if not (MIN_QUERY_LEN <= len(sequence) <= MAX_QUERY_LEN):
            return None
        seq = igem.validate_sequence(igem.clean_sequence(sequence))

        self.session.open_search_page()
        job_id = self.session.submit(seq)
        job = self.session.wait_job(job_id)

        data = job.get("data") or {}
        returned = igem.clean_sequence(str(data.get("sequence") or ""))
        query_length = len(returned) or len(seq)

        return_value = job.get("returnvalue") or {}
        if return_value.get("success") is False:
            raise LookupFailed(str(return_value.get("error") or "BLAST 返回失败。"))

        hits = [igem.Hit.from_api(item)
                for item in (return_value.get("hits") or [])
                if isinstance(item, dict)]
        best = _pick_exact(hits, query_length)
        return (best.part_name or None) if best is not None else None

    def query_many(self, sequences: Sequence[str],
                   progress: Optional[Callable[[str], None]] = None,
                   on_step: Optional[Callable[[int, int], None]] = None
                   ) -> Dict[str, Optional[str]]:
        """按序列去重后逐条查询，返回 {序列: 砖号或 None}。"""
        emit = progress if progress is not None else self._progress
        results: Dict[str, Optional[str]] = {}
        todo: List[str] = []

        for s in sequences:
            if not (MIN_QUERY_LEN <= len(s) <= MAX_QUERY_LEN):
                results.setdefault(s, None)
                continue
            if s in results or s in todo:
                continue
            if self.cache is not None:
                cached = self.cache.get(s)
                if cached is not SharedCache.MISS:
                    results[s] = cached
                    continue
            todo.append(s)

        if self.cache is not None and len(results):
            emit(f"有 {len(results)} 段直接用缓存结果")

        consecutive_errors = 0
        for index, seq in enumerate(todo, 1):
            if on_step is not None:
                try:
                    on_step(index, len(todo))
                except Exception:  # noqa: BLE001
                    pass
            emit(f"查砖号 {index}/{len(todo)}（{len(seq)} bp）…")
            brick: Optional[str] = None
            last_error: Optional[Exception] = None
            for attempt in (1, 2):
                try:
                    brick = self.query(seq)
                    last_error = None
                    break
                except Exception as exc:  # noqa: BLE001
                    last_error = exc
                    if attempt == 1:
                        emit(f"这一条出错了，重试一次：{exc}")
                        time.sleep(1.5)
            if last_error is not None:
                consecutive_errors += 1
                results[seq] = None
                if self.cache is not None:
                    self.cache.put(seq, None)
                emit(f"第 {index} 条查询失败，改用序列：{last_error}")
                if consecutive_errors >= 3:
                    raise LookupFailed(
                        "连续 3 条都查询失败，浏览器可能已经断开，已停止这个文件。\n"
                        f"最后一条的错误：{last_error}"
                    )
                continue
            consecutive_errors = 0
            results[seq] = brick
            if self.cache is not None:
                self.cache.put(seq, brick)

        for s in sequences:
            results.setdefault(s, None)
        return results


def preflight_login(profile_dir: Optional[Path] = None, *, headless: bool = True,
                    verbose: bool = False,
                    progress: Optional[Callable[[str], None]] = None) -> bool:
    """先单独确认一次登录态；过期就弹一次窗口让用户登录。返回是否新登录。"""
    base = Path(profile_dir) if profile_dir else igem.DEFAULT_PROFILE
    say = progress or (lambda _m: None)
    session = igem.Session(headless=headless, profile_dir=base, verbose=verbose)
    session.start()
    try:
        session.open_search_page()
        session.ensure_login(interactive=False)
        say("登录态有效。")
        return False
    except igem.LoginRequired:
        say("登录态已过期，打开浏览器窗口等你登录…")
    finally:
        session.close()

    session = igem.Session(headless=False, profile_dir=base, verbose=verbose)
    session.start()
    try:
        session.open_search_page()
        session.ensure_login(interactive=True)
        say("登录完成。")
        return True
    finally:
        session.close()


# ==========================================================================
# 4. 单文件流程
# ==========================================================================

def read_dna_file(path: Path) -> bytes:
    data = path.read_bytes()
    if data[:2] == b"\x1f\x8b":                 # gzip 压缩的 .dna
        data = gzip.decompress(data)
    return data


def default_output_path(dna_path: Path) -> Path:
    return dna_path.with_suffix(".txt")


class NoAnnotations(ValueError):
    """文件里一条标注都没有，没东西可以按段输出。"""


_SEGMENT_NAMES = {
    0x00: "DNA 序列", 0x02: "序列设置", 0x03: "酶列表", 0x05: "引物",
    0x06: "备注", 0x07: "序列锚点", 0x08: "序列属性", 0x09: "文件头",
    0x0A: "标注 Features", 0x0B: "序列历史", 0x0D: "自定义标签",
    0x0E: "自定义酶集", 0x10: "测序峰图", 0x11: "序列比对",
    0x13: "引物列表", 0x14: "引物位置", 0x15: "序列说明",
    0x1B: "比对数据 BAM", 0x1C: "酶显示设置",
}


def segment_inventory(segs: Sequence[Tuple[int, bytes]]) -> str:
    """把数据段清点成一行，方便看出文件里到底有什么。"""
    counts: Dict[int, int] = {}
    for seg_type, _payload in segs:
        counts[seg_type] = counts.get(seg_type, 0) + 1
    parts = []
    for seg_type in sorted(counts):
        name = _SEGMENT_NAMES.get(seg_type, f"未知(0x{seg_type:02x})")
        n = counts[seg_type]
        parts.append(f"{name}×{n}" if n > 1 else name)
    return "、".join(parts)


def describe_file(dna_path: Path, *, include_primers: bool = False) -> str:
    """把 .dna 里到底有什么列出来（排查用，不联网）。"""
    data = read_dna_file(dna_path)
    segs = read_segments(data)
    parsed = parse_snapgene(data, include_primers=include_primers)
    has_feature_seg = any(t == SEG_FEATURES for t, _ in segs)

    lines = [
        f"文件：{dna_path.name}（{len(data):,} 字节）",
        f"数据段：{len(segs)} 个 —— {segment_inventory(segs)}",
        f"序列：{len(parsed.seq):,} bp · "
        f"{'环状' if parsed.topology == 'circular' else '线状'}",
        f"标注：{parsed.feature_count} 条 · 引物：{parsed.primer_count} 个",
    ]
    if not has_feature_seg:
        lines.append("注意：文件里根本没有「标注 Features」数据段（0x0A），"
                     "所以一个部件都读不到。")
    for f in parsed.features[:40]:
        coords = ", ".join(f"{a}-{b}" for a, b in f.segs)
        lines.append(f"  · {f.label} [{coords}] "
                     f"{'反向' if f.strand < 0 else '正向'} ({f.type})")
    if len(parsed.features) > 40:
        lines.append(f"  … 还有 {len(parsed.features) - 40} 条")
    return "\n".join(lines)


def _no_features_message(dna_path: Path, parsed: Parsed) -> str:
    head = [f"{dna_path.name} 里一条标注都没有，没有东西可以按段输出。", ""]
    try:
        head.append(describe_file(dna_path))
    except Exception:  # noqa: BLE001
        head.append(f"序列：{len(parsed.seq):,} bp · "
                    f"{'环状' if parsed.topology == 'circular' else '线状'}")
    head.append("")
    if parsed.primer_count:
        head.append(
            f"文件里只有 {parsed.primer_count} 个引物结合位点，没有正式标注。"
            "引物不是砖，默认不当成段；真要看的话加 --primers。")
    else:
        head.append(
            "这类文件通常是测序/比对工程（只有序列和测序峰图），"
            "不是带部件标注的质粒图谱。\n"
            "如果 SnapGene 里能看到部件，把那个 .dna 发我吧——"
            "可能是保存时没把标注带上。")
    head.append("确实只想把整条序列当一段输出，加 --force。")
    return "\n".join(head)


def plan(dna_path: Path, *, include_primers: bool = False,
         rev_comp: bool = True, require_features: bool = True,
         log: Optional[Callable[[str], None]] = None) -> Tuple[Parsed, List[Piece]]:
    """解析 .dna 并分段（不联网）。"""
    say = log or (lambda _m: None)
    parsed = parse_snapgene(read_dna_file(dna_path), include_primers=include_primers)
    say(f"读到质粒：{len(parsed.seq)} bp · "
        f"{'环状' if parsed.topology == 'circular' else '线状'} · "
        f"{parsed.feature_count} 条标注"
        + (f"（另有 {parsed.primer_count} 个引物）" if parsed.primer_count else ""))

    # 一条标注都没有：直接说清楚，别输出一整条"未标记序列"糊弄人
    if require_features and not parsed.features:
        raise NoAnnotations(_no_features_message(dna_path, parsed))

    pieces = build_pieces(parsed, rev_comp=rev_comp, log=say)
    if not pieces:
        raise ValueError("这条质粒上没有任何标注，也没有序列可输出。")

    covered_bases = sum(_ranges_len(p.ranges) for p in pieces)
    if covered_bases != len(parsed.seq):
        raise ValueError(
            f"内部错误：只覆盖了 {covered_bases}/{len(parsed.seq)} bp，会漏序列，已中止。")
    say(f"分成 {len(pieces)} 段，覆盖 {covered_bases}/{len(parsed.seq)} bp（完整覆盖）")
    return parsed, pieces


def apply_bricks(pieces: Sequence[Piece], found: Dict[str, Optional[str]],
                 log: Optional[Callable[[str], None]] = None) -> None:
    say = log or (lambda _m: None)
    total = len(pieces)
    for i, p in enumerate(pieces, 1):
        p.brick = found.get(p.seq)
        say(f"[{i}/{total}] {p.title} → {p.brick or '查不到砖号，写序列'}")


def write_txt(dna_path: Path, pieces: Sequence[Piece],
              out_path: Optional[Path] = None) -> Path:
    target = Path(out_path) if out_path else default_output_path(dna_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lines = output_lines(pieces)
    with open(target, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return target


def convert_file(dna_path: Path, *, finder: BrickFinder,
                 out_path: Optional[Path] = None,
                 include_primers: bool = False, rev_comp: bool = True,
                 require_features: bool = True,
                 log: Optional[Callable[[str], None]] = None,
                 on_step: Optional[Callable[[int, int], None]] = None
                 ) -> Tuple[Path, int]:
    """用已有的浏览器会话把一个 .dna 转成 txt。"""
    say = log or (lambda _m: None)
    _parsed, pieces = plan(dna_path, include_primers=include_primers,
                           rev_comp=rev_comp, require_features=require_features,
                           log=say)
    found = finder.query_many([p.seq for p in pieces], progress=say,
                              on_step=on_step)
    apply_bricks(pieces, found, log=say)
    target = write_txt(dna_path, pieces, out_path)
    say(f"已写出 {target.name}（{len(pieces)} 行）")
    return target, len(pieces)


def convert(dna_path: Path, *, out_path: Optional[Path] = None,
            include_primers: bool = False, rev_comp: bool = True,
            headless: bool = True, verbose: bool = False,
            timeout: float = 180.0, profile_dir: Optional[Path] = None,
            log: Optional[Callable[[str], None]] = None,
            dry_run: bool = False, force: bool = False) -> Tuple[Path, int]:
    """一个文件的完整流程（自己开浏览器）。"""
    require = not force
    if dry_run:
        _parsed, pieces = plan(dna_path, include_primers=include_primers,
                               rev_comp=rev_comp, require_features=require, log=log)
        return (Path(out_path) if out_path else default_output_path(dna_path)), len(pieces)
    with BrickFinder(profile_dir=profile_dir, headless=headless, verbose=verbose,
                     timeout=timeout, progress=log) as finder:
        return convert_file(dna_path, finder=finder, out_path=out_path,
                            include_primers=include_primers, rev_comp=rev_comp,
                            require_features=require, log=log)


# ==========================================================================
# 5. 随时丢文件进来的查询池
# ==========================================================================

def _worker_profile(base: Path, index: int) -> Path:
    """第 0 个用原来的配置目录，其余用兄弟目录。

    放在同一层是为了让 Session 去 profile_dir.parent/login-state.json
    找登录态时，仍然指向 app/login-state.json。
    """
    if index == 0:
        return base
    return base.parent / f"{base.name}_w{index + 1}"


IDLE_RELEASE_SECONDS = 90.0      # 闲这么久就把浏览器关掉，别白占内存


class Pool:
    """常驻的查询池：随时把文件丢进来，随时开查。

    · 有空位就立刻开查，没有空位就在队列里等着；
    · 每个 worker 自己一个浏览器（各自一个配置目录），只在有活时才开，
      闲下来一阵子自动关掉；
    · 同一个砖在几个文件里出现时只查一次。
    """

    def __init__(self, jobs: int = DEFAULT_JOBS, on_event: Optional[Event] = None, *,
                 out_dir: Optional[Path] = None,
                 include_primers: bool = False, rev_comp: bool = True,
                 headless: bool = True, timeout: float = 180.0,
                 profile_dir: Optional[Path] = None,
                 verbose: bool = False, force: bool = False) -> None:
        self.jobs = max(1, int(jobs or 1))
        self.emit: Event = on_event or (lambda *_a, **_k: None)
        self.out_dir = Path(out_dir) if out_dir else None
        self.include_primers = include_primers
        self.rev_comp = rev_comp
        self.force = force
        self.headless = headless
        self.timeout = timeout
        self.verbose = verbose
        self.base = Path(profile_dir) if profile_dir else igem.DEFAULT_PROFILE
        self.cache = SharedCache()

        self._todo: "queue.Queue" = queue.Queue()
        self._threads: List[threading.Thread] = []
        self._lock = threading.Lock()
        self._idle = threading.Condition(self._lock)
        self._ready_lock = threading.Lock()
        self._results: Dict[Path, Tuple[str, str]] = {}
        self._active = 0
        self._ready = False
        self._closed = False

    # -- 对外 -------------------------------------------------------------

    def start(self) -> None:
        with self._lock:
            if self._threads or self._closed:
                return
            for index in range(self.jobs):
                t = threading.Thread(target=self._worker, args=(index,), daemon=True)
                self._threads.append(t)
                t.start()

    def submit(self, path: Path, out_path: Optional[Path] = None) -> None:
        """丢一个文件进来：有空位就现在开查，没空位就排队。"""
        self.start()
        with self._lock:
            if self._closed:
                raise LookupFailed("查询池已经关闭了。")
            self._active += 1
        self._todo.put((Path(path), Path(out_path) if out_path else None))

    def wait_idle(self) -> Dict[Path, Tuple[str, str]]:
        """等所有已经丢进来的文件都查完。"""
        with self._idle:
            while self._active > 0:
                self._idle.wait(1.0)
            return dict(self._results)

    def shutdown(self, timeout: float = 6.0) -> None:
        with self._lock:
            if self._closed:
                return
            self._closed = True
            threads = list(self._threads)
        for _ in threads:
            self._todo.put(None)
        deadline = time.time() + timeout
        for t in threads:
            t.join(timeout=max(0.1, deadline - time.time()))
        for index in range(1, self.jobs):
            profile = _worker_profile(self.base, index)
            if profile.exists():
                shutil.rmtree(profile, ignore_errors=True)

    # -- 内部 -------------------------------------------------------------

    def _ensure_ready(self) -> None:
        """整个池只确认一次登录态：同时开多个浏览器时这是必须的。"""
        with self._ready_lock:
            if self._ready:
                return
            if self.jobs > 1:
                preflight_login(self.base, headless=self.headless,
                                verbose=self.verbose,
                                progress=lambda m: self.emit("log", None, m))
                # 让刚关掉的浏览器把配置目录的锁松开
                time.sleep(0.8)
            self._ready = True

    def _finish_one(self) -> None:
        with self._idle:
            self._active -= 1
            self._idle.notify_all()

    def _worker(self, index: int) -> None:
        finder: Optional[BrickFinder] = None
        try:
            while True:
                try:
                    item = self._todo.get(timeout=IDLE_RELEASE_SECONDS)
                except queue.Empty:
                    # 闲了一阵子：把浏览器关掉，别白占内存
                    if finder is not None:
                        finder.close()
                        finder = None
                    continue
                if item is None:
                    return
                path, out_path = item
                self.emit("start", path, None)
                try:
                    self._ensure_ready()
                    if finder is None:
                        finder = BrickFinder(
                            profile_dir=_worker_profile(self.base, index),
                            headless=self.headless, verbose=self.verbose,
                            timeout=self.timeout, cache=self.cache,
                            allow_interactive_login=(self.jobs == 1),
                        )
                        finder.open()
                    finder.set_progress(lambda m, _p=path: self.emit("log", _p, m))
                    target, count = convert_file(
                        path, finder=finder,
                        out_path=((self.out_dir / f"{path.stem}.txt")
                                  if self.out_dir else out_path),
                        include_primers=self.include_primers,
                        rev_comp=self.rev_comp,
                        require_features=not self.force,
                        log=lambda m, _p=path: self.emit("log", _p, m),
                        on_step=lambda i, n, _p=path: self.emit("step", _p, (i, n)),
                    )
                    with self._lock:
                        self._results[path] = ("done", f"{target}（{count} 行）")
                    self.emit("done", path, (target, count))
                except Exception as exc:  # noqa: BLE001
                    with self._lock:
                        self._results[path] = ("fail", str(exc))
                    self.emit("fail", path, exc)
                    # 浏览器可能已经挂了，下一个文件重新开一个
                    if finder is not None:
                        finder.close()
                        finder = None
                finally:
                    self._finish_one()
        finally:
            if finder is not None:
                finder.close()


def run_files(paths: Sequence[Path], *, jobs: int = DEFAULT_JOBS,
              on_event: Optional[Event] = None, **opts) -> Dict[Path, Tuple[str, str]]:
    """命令行用：把所有文件丢进池子，等它们全部查完。"""
    paths = [Path(p) for p in paths]
    # 只有一个文件时没必要开多个浏览器，也没必要先跑一遍登录确认
    pool = Pool(jobs=max(1, min(int(jobs or 1), len(paths) or 1)),
                on_event=on_event, **opts)
    try:
        for p in paths:
            pool.submit(p)
        return pool.wait_idle()
    finally:
        pool.shutdown()


# ==========================================================================
# 6. 界面
# ==========================================================================

# 配色：板岩 / 紫 / 蓝 / 青 四组，每组深—中—浅三档
INK       = "#465E77"   # 主文字、深色按钮
INK_SOFT  = "#7C8B9B"   # 次要文字
GRAY      = "#C1C8CE"   # 中性
MIST      = "#E1E7EC"   # 底色、分隔线
PAPER     = "#FFFFFF"
CANVAS    = "#F6F9FB"   # 日志底
VIOLET    = "#8F8BD9"
VIOLET_L  = "#C7C6F2"
VIOLET_XL = "#E0E3FB"
BLUE      = "#4C8ACA"
BLUE_L    = "#95C3F1"
BLUE_XL   = "#D2E6FB"
TEAL      = "#1D9695"
TEAL_L    = "#9CDCDD"
TEAL_XL   = "#D5F3F4"
ROSE      = "#BE5A6A"   # 只用在报错
ROSE_XL   = "#F7E7EA"

FONT_CN    = "Microsoft YaHei UI"
FONT_UI    = (FONT_CN, 9)
FONT_UI_B  = (FONT_CN, 9, "bold")
FONT_TINY  = (FONT_CN, 8)
FONT_TITLE = (FONT_CN, 16, "bold")
FONT_NAME  = (FONT_CN, 11, "bold")
FONT_MONO  = ("Consolas", 9)

GUI_MAX_PARALLEL = 3
ICON_PATH = HERE / "fish.ico"


def run_gui() -> None:
    import tkinter as tk
    from tkinter import filedialog, ttk

    def flat_button(parent, text: str, command, *, bg: str, fg: str, hover: str,
                    font=FONT_UI_B, padx: int = 16, pady: int = 8):
        """扁平按钮：没有立体边框，鼠标移上去换底色。"""
        btn = tk.Button(parent, text=text, command=command, font=font,
                        bg=bg, fg=fg, activebackground=hover, activeforeground=fg,
                        bd=0, relief="flat", highlightthickness=0,
                        cursor="hand2", padx=padx, pady=pady)

        def on_enter(_e=None):
            if str(btn["state"]) != "disabled":
                btn.configure(bg=hover)

        def on_leave(_e=None):
            if str(btn["state"]) != "disabled":
                btn.configure(bg=bg)

        def recolor(new_bg: str, new_fg: str, new_hover: str) -> None:
            nonlocal bg, fg, hover
            bg, fg, hover = new_bg, new_fg, new_hover
            btn.configure(bg=bg, fg=fg, activebackground=hover, activeforeground=fg)

        btn.bind("<Enter>", on_enter)
        btn.bind("<Leave>", on_leave)
        btn.recolor = recolor           # type: ignore[attr-defined]
        return btn

    class ThinBar(tk.Frame):
        """4px 高的进度条。"""

        def __init__(self, parent, *, color: str = BLUE):
            super().__init__(parent, bg=MIST, height=4, bd=0, highlightthickness=0)
            self.pack_propagate(False)
            self._fill = tk.Frame(self, bg=color, bd=0, highlightthickness=0)
            self._fill.place(x=0, y=0, relwidth=0.0, relheight=1.0)

        def set(self, fraction: float) -> None:
            self._fill.place_configure(relwidth=max(0.0, min(1.0, float(fraction))))

        def color(self, value: str) -> None:
            self._fill.configure(bg=value)

    root = tk.Tk()
    root.title("质粒 → 砖号")
    root.configure(bg=PAPER)
    try:
        if ICON_PATH.exists():
            root.iconbitmap(default=str(ICON_PATH))     # 一笔画小鱼
    except Exception:  # noqa: BLE001
        pass

    style = ttk.Style()
    try:
        style.theme_use("clam")
    except tk.TclError:
        pass
    style.configure("Slim.Vertical.TScrollbar", gripcount=0, relief="flat",
                    background=GRAY, darkcolor=PAPER, lightcolor=PAPER,
                    troughcolor=PAPER, bordercolor=PAPER, arrowcolor=INK_SOFT,
                    width=9)
    style.map("Slim.Vertical.TScrollbar", background=[("active", INK_SOFT)])

    # ---------------- 主窗口 ----------------

    header = tk.Frame(root, bg=VIOLET, height=64)
    header.pack(fill="x")
    header.pack_propagate(False)
    tk.Label(header, text="质粒 → 砖号", bg=VIOLET, fg="#FFFFFF",
             font=FONT_TITLE).pack(side="left", padx=26)

    wrap = tk.Frame(root, bg=PAPER, padx=26, pady=22)
    wrap.pack(fill="both", expand=True)

    tk.Label(wrap, text="输入 SnapGene .dna", bg=PAPER, fg=INK,
             font=(FONT_CN, 10), anchor="w").pack(fill="x")

    pick_btn = flat_button(wrap, "选择质粒文件…", lambda: on_pick(),
                           bg=BLUE, fg="#FFFFFF", hover=INK, pady=10)
    pick_btn.pack(fill="x", pady=(10, 0))

    tk.Label(wrap, text="可同时查找多个文件", bg=PAPER, fg=INK_SOFT,
             font=FONT_TINY, anchor="w").pack(fill="x", pady=(14, 0))

    events: "queue.Queue" = queue.Queue()
    windows: Dict[Path, "FileWindow"] = {}
    holder: Dict[str, Optional[Pool]] = {"pool": None}

    def on_event(kind: str, path, payload) -> None:
        events.put((kind, path, payload))

    def get_pool() -> Pool:
        if holder["pool"] is None:
            holder["pool"] = Pool(GUI_MAX_PARALLEL, on_event=on_event)
        return holder["pool"]

    # ---------------- 每个质粒一个窗口 ----------------

    class FileWindow:
        def __init__(self, path: Path) -> None:
            self.path = path
            self.kind = "queued"       # queued / running / done / failed
            self.target: Optional[Path] = None

            self.top = tk.Toplevel(root)
            self.top.title(path.name)
            self.top.geometry("760x520")
            self.top.minsize(560, 380)
            self.top.configure(bg=PAPER)

            band = tk.Frame(self.top, bg=TEAL)
            band.pack(fill="x")
            inner = tk.Frame(band, bg=TEAL, padx=24, pady=12)
            inner.pack(fill="both", expand=True)
            tk.Label(inner, text=path.name, bg=TEAL, fg="#FFFFFF", font=FONT_NAME,
                     anchor="w", justify="left", wraplength=690).pack(fill="x")
            tk.Label(inner, text=str(path.parent), bg=TEAL, fg=TEAL_XL,
                     font=FONT_TINY, anchor="w").pack(fill="x", pady=(3, 0))

            row = tk.Frame(self.top, bg=PAPER, padx=24)
            row.pack(fill="x", pady=(14, 8))
            self.pill = tk.Label(row, text="排队中", bg=MIST, fg=INK,
                                 font=FONT_UI_B, padx=12, pady=4)
            self.pill.pack(side="left")
            self.step_label = tk.Label(row, text="等前面的查完", bg=PAPER,
                                       fg=INK_SOFT, font=FONT_TINY)
            self.step_label.pack(side="right")

            bar_holder = tk.Frame(self.top, bg=PAPER, padx=24)
            bar_holder.pack(fill="x")
            self.bar = ThinBar(bar_holder, color=BLUE)
            self.bar.pack(fill="x", pady=(4, 14))

            # 先把底部按钮占住，再让日志区吃掉剩下的空间（顺序反了按钮会被挤出去）
            foot = tk.Frame(self.top, bg=PAPER, padx=24, pady=14)
            foot.pack(side="bottom", fill="x")
            self.open_btn = flat_button(foot, "打开 txt", self._open_txt,
                                        bg=MIST, fg=INK, hover=BLUE_XL,
                                        font=FONT_UI, padx=14, pady=7)
            self.open_btn.configure(state="disabled", disabledforeground=GRAY)
            self.open_btn.pack(side="right")

            body = tk.Frame(self.top, bg=PAPER, padx=24)
            body.pack(fill="both", expand=True)
            self.log_box = tk.Text(body, font=FONT_MONO, bg=CANVAS, fg="#4A5968",
                                   bd=0, highlightthickness=1, highlightbackground=MIST,
                                   highlightcolor=MIST, padx=12, pady=10,
                                   wrap="word", state="disabled")
            sb = ttk.Scrollbar(body, orient="vertical", command=self.log_box.yview,
                               style="Slim.Vertical.TScrollbar")
            self.log_box.configure(yscrollcommand=sb.set)
            sb.pack(side="right", fill="y")
            self.log_box.pack(side="left", fill="both", expand=True)

        # -- 状态 ---------------------------------------------------------

        def append(self, msg: str) -> None:
            self.log_box.configure(state="normal")
            self.log_box.insert("end", msg + "\n")
            self.log_box.see("end")
            self.log_box.configure(state="disabled")

        def set_queued(self) -> None:
            self.kind = "queued"
            self.pill.configure(text="排队中", bg=MIST, fg=INK)
            self.step_label.configure(text="等前面的查完")
            self.bar.set(0.0)
            self.top.title(self.path.name)

        def set_running(self) -> None:
            self.kind = "running"
            self.pill.configure(text="查询中", bg=BLUE_XL, fg=BLUE)
            self.step_label.configure(text="准备查询…")
            self.top.title(self.path.name)

        def set_step(self, index: int, total: int) -> None:
            if total <= 0:
                return
            self.step_label.configure(text=f"第 {index}/{total} 段")
            self.bar.set(index / total)

        def finish(self, target: Path, count: int = 0) -> None:
            self.kind = "done"
            self.target = target
            self.pill.configure(text="完成", bg=TEAL_XL, fg=TEAL)
            self.step_label.configure(text=(f"{count} 行 · 已生成 txt" if count
                                            else "已生成 txt"))
            self.bar.set(1.0)
            self.top.title("✓ " + self.path.name)
            self.open_btn.recolor("#1D9695", "#FFFFFF", INK)   # type: ignore[attr-defined]
            self.open_btn.configure(state="normal")

        def fail(self, exc: Exception) -> None:
            self.kind = "failed"
            self.pill.configure(text="失败", bg=ROSE_XL, fg=ROSE)
            self.top.title("✗ " + self.path.name)

        def _open_txt(self) -> None:
            if self.target and self.target.exists():
                try:
                    os.startfile(str(self.target))  # type: ignore[attr-defined]
                except Exception as exc:  # noqa: BLE001
                    self.append(f"打不开：{exc}")

    # ---------------- 交互 ----------------

    def on_pick() -> None:
        names = filedialog.askopenfilenames(
            title="选择质粒文件",
            filetypes=[("SnapGene 质粒文件", "*.dna *.snapgene"), ("所有文件", "*.*")],
        )
        for name in names:
            path = Path(name).resolve()
            existing = windows.get(path)
            if existing is not None:
                existing.top.lift()
                continue
            win = FileWindow(path)
            windows[path] = win
            try:
                get_pool().submit(path)
            except Exception as exc:  # noqa: BLE001
                win.fail(exc)

    def drain() -> None:
        try:
            while True:
                kind, path, payload = events.get_nowait()
                if kind == "log":
                    if path is None:
                        for win in windows.values():
                            if win.kind in ("queued", "running"):
                                win.append(payload)
                    else:
                        win = windows.get(path)
                        if win is not None:
                            win.append(payload)
                elif kind == "step":
                    win = windows.get(path)
                    if win is not None:
                        win.set_step(*payload)
                elif kind == "start":
                    win = windows.get(path)
                    if win is not None:
                        win.set_running()
                elif kind == "done":
                    target, count = payload
                    win = windows.get(path)
                    if win is not None:
                        win.append("")
                        win.append(f"共 {count} 行 → {target}")
                        win.finish(target, count)
                elif kind == "fail":
                    win = windows.get(path)
                    if win is not None:
                        win.append("")
                        win.append(str(payload))
                        win.fail(payload)
        except queue.Empty:
            pass
        root.after(120, drain)

    def on_close() -> None:
        if holder["pool"] is not None:
            holder["pool"].shutdown()
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)
    root.update_idletasks()
    root.geometry(f"{max(520, root.winfo_reqwidth())}x{root.winfo_reqheight()}")
    root.resizable(False, False)
    root.after(120, drain)
    root.mainloop()


# ==========================================================================
# 7. 命令行
# ==========================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dna_to_bricks",
        description="把质粒 .dna 转成 txt：一行一段，「名字：砖号」；"
                    "查不到砖号的那一行写「名字：序列」。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            '  python dna_to_bricks.py 质粒.dna                # 生成 质粒.txt\n'
            '  python dna_to_bricks.py a.dna b.dna c.dna       # 三个一起查\n'
            '  python dna_to_bricks.py *.dna -j 4              # 同时查 4 个\n'
            "  python dna_to_bricks.py --gui                   # 打开窗口\n"
            "  python dna_to_bricks.py 质粒.dna --dry-run      # 只看分段，不联网\n"
        ),
    )
    parser.add_argument("dna", nargs="*", help="质粒文件路径（.dna / .snapgene），可以给多个")
    parser.add_argument("-o", "--out", help="输出 txt 路径（只在单文件时用；默认同名 .txt）")
    parser.add_argument("-j", "--jobs", type=int, default=DEFAULT_JOBS,
                        help=f"同时处理几个文件，默认 {DEFAULT_JOBS}")
    parser.add_argument("--gui", action="store_true", help="打开窗口")
    parser.add_argument("--primers", action="store_true",
                        help="把引物结合位点也当成一段（默认不算，引物不是砖）")
    parser.add_argument("--no-revcomp", action="store_true",
                        help="反向链的段不取反向互补，直接按正链输出")
    parser.add_argument("--show", action="store_true",
                        help="显示浏览器窗口（默认后台跑，需要登录时会自己弹出来）")
    parser.add_argument("--timeout", type=float, default=180.0,
                        help="等待单条 BLAST 结果的秒数，默认 180")
    parser.add_argument("--profile", help=f"登录态目录（默认 {igem.DEFAULT_PROFILE}）")
    parser.add_argument("--dry-run", action="store_true",
                        help="只解析分段并打印，不联网查砖号")
    parser.add_argument("--info", action="store_true",
                        help="只打印这个 .dna 里到底有什么（数据段/标注/坐标），不联网")
    parser.add_argument("--force", action="store_true",
                        help="文件里没有标注时也照样输出（整条序列当一段）")
    parser.add_argument("-v", "--verbose", action="store_true", help="打印过程信息")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    if sys.stdout is None:
        sys.stdout = open(os.devnull, "w", encoding="utf-8")
    if sys.stderr is None:
        sys.stderr = open(os.devnull, "w", encoding="utf-8")
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass

    args = build_parser().parse_args(argv)

    if args.gui or not args.dna:
        run_gui()
        return 0

    paths = [Path(p) for p in args.dna]
    missing = [p for p in paths if not p.exists()]
    if missing:
        for p in missing:
            print(f"找不到文件：{p}", file=sys.stderr)
        return 1

    def say(msg: str) -> None:
        print(msg, file=sys.stderr, flush=True)

    common = dict(include_primers=args.primers, rev_comp=not args.no_revcomp,
                  headless=not args.show, verbose=args.verbose,
                  timeout=args.timeout, force=args.force,
                  profile_dir=Path(args.profile) if args.profile else None)

    # --info：只把这个文件里有什么列出来
    if args.info:
        for path in paths:
            try:
                print(describe_file(path, include_primers=args.primers))
            except Exception as exc:  # noqa: BLE001
                print(f"{path.name}: 读不出来 —— {exc}")
            print()
        return 0

    # --dry-run：只分段
    if args.dry_run:
        for path in paths:
            try:
                _parsed, pieces = plan(path, include_primers=common["include_primers"],
                                       rev_comp=common["rev_comp"],
                                       require_features=not args.force, log=say)
            except Exception as exc:  # noqa: BLE001
                print(f"失败：{exc}", file=sys.stderr)
                continue
            print(f"# {path.name}: {len(pieces)} 段（dry-run，未联网）", file=sys.stderr)
            for i, p in enumerate(pieces, 1):
                print(f"{i:>3}. {p.title}  {p.coords}  {len(p.seq)} bp", file=sys.stderr)
        return 0

    # 单文件且指定了 -o：走单文件流程
    if len(paths) == 1 and args.out:
        try:
            target, count = convert(paths[0], out_path=Path(args.out), log=say, **common)
        except KeyboardInterrupt:
            print("已中断。", file=sys.stderr)
            return 130
        except Exception as exc:  # noqa: BLE001
            print(f"失败：{exc}", file=sys.stderr)
            if args.verbose:
                traceback.print_exc()
            return 1
        print(target)
        return 0

    def on_event(kind: str, path, payload) -> None:
        if kind == "log":
            print(("[总] " if path is None else f"[{path.name}] ") + str(payload),
                  file=sys.stderr, flush=True)
        elif kind == "start":
            print(f"[{path.name}] 开始", file=sys.stderr, flush=True)
        elif kind == "done":
            target, count = payload
            print(f"[{path.name}] 完成 {target}（{count} 行）",
                  file=sys.stderr, flush=True)
        elif kind == "fail":
            print(f"[{path.name}] 失败：{payload}", file=sys.stderr, flush=True)

    try:
        results = run_files(paths, jobs=args.jobs, on_event=on_event, **common)
    except KeyboardInterrupt:
        print("已中断。", file=sys.stderr)
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"失败：{exc}", file=sys.stderr)
        if args.verbose:
            traceback.print_exc()
        return 1

    ok = sum(1 for s, _ in results.values() if s == "done")
    for path in paths:
        status, detail = results.get(path, ("fail", "没有结果"))
        if status == "done":
            print(detail.split("（")[0])
    print(f"完成 {ok}/{len(paths)} 个。", file=sys.stderr)
    return 0 if ok == len(paths) else 1


if __name__ == "__main__":
    raise SystemExit(main())
