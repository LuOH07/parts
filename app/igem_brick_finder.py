#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
iGEM 砖号查询器 (igem_brick_finder)
====================================

给定一段 DNA 序列，自动完成下面这套人工操作：

  1. 打开 https://registry.igem.org/parts/blast
  2. 在 "Sequence" 文本框中填入序列，点击 "Search"
  3. 等待页面自动跳转到结果页 /parts/blast/<jobId>
  4. 点击结果页上的 "Full match" 过滤按钮
  5. 只看过滤后的第一个结果：
        - 相似度 100%
        - 且长度与输入序列完全一致
     满足则输出该砖的砖号（如 BBa_E0040）

依赖：playwright（首次运行需要联网安装浏览器，见 README.md）
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import threading
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional
from urllib.parse import quote

# --------------------------------------------------------------------------
# 站点常量
# --------------------------------------------------------------------------

SITE_URL = "https://registry.igem.org/parts/blast"
API_BASE = "https://api.registry.igem.org/v1"
ME_URL = API_BASE + "/accounts/me"
JOB_URL = API_BASE + "/blast/jobs/{job_id}"

APP_VERSION = "2.2"

ALGORITHMS = ("megablast", "blastn")
DEFAULT_ALGORITHM = "megablast"

NUCLEOTIDE_RE = re.compile(r"^[ACGTURYKMSWBDHVNX-]+$", re.IGNORECASE)
JOB_ID_RE = re.compile(r"/parts/blast/([0-9A-Za-z_.:-]+)$")

# 站点自身 "Full match" 的判定阈值（见前端源码）
FULL_QUERY_COVERAGE = 99.5
IDENTITY_EXACT = 99.999

if getattr(sys, "frozen", False):
    # PyInstaller 打包后，__file__ 在临时目录里；配置要放在 exe 旁边
    HERE = Path(sys.executable).resolve().parent
else:
    HERE = Path(__file__).resolve().parent
DEFAULT_PROFILE = HERE / ".igem_profile"


# --------------------------------------------------------------------------
# 数据模型
# --------------------------------------------------------------------------

@dataclass
class Hit:
    """BLAST 命中（字段名与 registry 前端接口一致）。"""
    part_name: str
    identity: float
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    subject_length: int
    alignment_length: int = 0
    mismatches: int = 0
    gaps: int = 0
    evalue: float = 0.0
    bit_score: float = 0.0

    @classmethod
    def from_api(cls, raw: Dict[str, Any]) -> "Hit":
        def as_int(key: str, default: int = 0) -> int:
            try:
                return int(raw.get(key) or default)
            except (TypeError, ValueError):
                return default

        def as_float(key: str, default: float = 0.0) -> float:
            try:
                return float(raw.get(key) if raw.get(key) is not None else default)
            except (TypeError, ValueError):
                return default

        return cls(
            part_name=str(raw.get("partName") or "").strip(),
            identity=as_float("identity"),
            query_start=as_int("queryStart", 1),
            query_end=as_int("queryEnd", 0),
            subject_start=as_int("subjectStart", 1),
            subject_end=as_int("subjectEnd", 0),
            subject_length=as_int("subjectLength", 0),
            alignment_length=as_int("alignmentLength", 0),
            mismatches=as_int("mismatches", 0),
            gaps=as_int("gaps", 0),
            evalue=as_float("evalue"),
            bit_score=as_float("bitScore"),
        )


@dataclass
class Result:
    status: str                      # exact / not_exact / no_hits / error
    part: Optional[str] = None
    query_length: int = 0
    identity: Optional[float] = None
    subject_length: Optional[int] = None
    job_id: Optional[str] = None
    url: Optional[str] = None
    message: str = ""
    first_hit: Optional[Hit] = None
    hits_total: int = 0
    hits_full_match: int = 0
    source: str = "api"              # api / dom

    def to_dict(self) -> Dict[str, Any]:
        data = {
            "status": self.status,
            "part": self.part,
            "query_length": self.query_length,
            "identity": self.identity,
            "subject_length": self.subject_length,
            "job_id": self.job_id,
            "url": self.url,
            "message": self.message,
            "hits_total": self.hits_total,
            "hits_full_match": self.hits_full_match,
            "source": self.source,
        }
        if self.first_hit is not None:
            data["first_hit"] = {
                "part": self.first_hit.part_name,
                "identity": self.first_hit.identity,
                "query_start": self.first_hit.query_start,
                "query_end": self.first_hit.query_end,
                "subject_start": self.first_hit.subject_start,
                "subject_end": self.first_hit.subject_end,
                "subject_length": self.first_hit.subject_length,
                "alignment_length": self.first_hit.alignment_length,
                "mismatches": self.first_hit.mismatches,
                "gaps": self.first_hit.gaps,
                "evalue": self.first_hit.evalue,
                "bit_score": self.first_hit.bit_score,
            }
        return data


class LoginRequired(RuntimeError):
    """没有登录 iGEM，无法使用 BLAST。"""


class SearchFailed(RuntimeError):
    """搜索过程出错。"""


# --------------------------------------------------------------------------
# 序列处理
# --------------------------------------------------------------------------

def clean_sequence(raw: str) -> str:
    """按站点的方式清洗序列：去掉 FASTA 头、空白，并转成大写。"""
    text = raw.replace("\r\n", "\n").replace("\r", "\n")
    text = re.sub(r"^>.*\n?", "", text, count=1)
    text = re.sub(r"\s+", "", text)
    return text.upper()


def validate_sequence(seq: str) -> str:
    if not seq:
        raise ValueError("序列为空。")
    if len(seq) < 7:
        raise ValueError(f"序列太短（{len(seq)} bp），至少需要 7 bp。")
    if len(seq) > 100000:
        raise ValueError(f"序列太长（{len(seq)} bp），上限 100000 bp。")
    if not NUCLEOTIDE_RE.match(seq):
        bad = sorted({c for c in seq if not NUCLEOTIDE_RE.match(c)})
        raise ValueError(f"序列里有非核苷酸字符：{''.join(bad)!r}"
                         "（只支持 A/C/G/T/U 及 IUPAC 简并碱基）。")
    return seq


# --------------------------------------------------------------------------
# 命中判定（与站点前端逻辑保持一致）
# --------------------------------------------------------------------------

def _bounds(hit: Hit):
    lo = min(hit.subject_start, hit.subject_end)
    hi = max(hit.subject_start, hit.subject_end)
    return lo, hi


def is_full_match(hit: Hit, query_length: int) -> bool:
    """站点 "Full match" 按钮的筛选条件。"""
    lo, hi = _bounds(hit)
    full_subject = lo == 1 and hi == hit.subject_length
    if query_length <= 0:
        return False
    query_coverage = (hit.query_end - hit.query_start + 1) / query_length * 100.0
    return full_subject and query_coverage >= FULL_QUERY_COVERAGE


def is_strict_exact(hit: Hit, query_length: int) -> bool:
    """用户要求：相似度 100% 且长度完全一致。"""
    return (
        hit.identity >= IDENTITY_EXACT
        and hit.subject_length == query_length
        and hit.query_start == 1
        and hit.query_end == query_length
        and hit.gaps == 0
        and hit.mismatches == 0
    )


def judge(hits: List[Hit], query_length: int) -> Dict[str, Any]:
    """在 "Full match" 结果中取第一条，判断是不是完全一致的砖。"""
    full = [h for h in hits if is_full_match(h, query_length)]
    total = len(hits)
    if not full:
        return {
            "status": "no_hits" if total == 0 else "no_full_match",
            "first": None,
            "full_count": 0,
            "total": total,
            "message": "没有命中任何部件。" if total == 0
                       else "没有任何结果能通过 Full match 过滤。",
        }

    first = full[0]
    if is_strict_exact(first, query_length):
        return {
            "status": "exact",
            "first": first,
            "full_count": len(full),
            "total": total,
            "message": f"第一个 Full match 结果是 {first.part_name}，"
                       f"相似度 100%，长度 {first.subject_length} bp 与输入一致。",
        }
    return {
        "status": "not_exact",
        "first": first,
        "full_count": len(full),
        "total": total,
        "message": (f"第一个 Full match 结果是 {first.part_name}，"
                    f"但相似度 {first.identity:.2f}%、长度 {first.subject_length} bp"
                    f"（输入 {query_length} bp），并非 100% 等长匹配。"),
    }


# --------------------------------------------------------------------------
# 浏览器会话
# --------------------------------------------------------------------------

def _import_playwright():
    try:
        from playwright.sync_api import sync_playwright  # noqa: F401
        return sync_playwright
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "缺少 playwright。请先安装：\n"
            "  python -m pip install playwright\n"
            "  python -m playwright install chromium"
        ) from exc


class Session:
    """一次带持久化登录态的浏览器会话。"""

    def __init__(
        self,
        *,
        headless: bool = True,
        profile_dir: Optional[Path] = None,
        algorithm: str = DEFAULT_ALGORITHM,
        timeout: float = 180.0,
        login_timeout: float = 300.0,
        verbose: bool = False,
        debug_dir: Optional[Path] = None,
        page_router: Optional[Callable[[Any], None]] = None,
        assume_logged_in: bool = False,
        progress: Optional[Callable[[str], None]] = None,
    ) -> None:
        self.headless = headless
        self.profile_dir = Path(profile_dir or DEFAULT_PROFILE)
        self.algorithm = algorithm if algorithm in ALGORITHMS else DEFAULT_ALGORITHM
        self.timeout_ms = int(timeout * 1000)
        self.login_timeout = login_timeout
        self.verbose = verbose
        self.debug_dir = debug_dir
        self.page_router = page_router
        self.assume_logged_in = assume_logged_in
        self.progress = progress

        self._pw = None
        self.ctx = None
        self.page = None
        self._job_responses: List[Any] = []
        self._watched_pages: set = set()

    # -- 生命周期 ---------------------------------------------------------

    def __enter__(self) -> "Session":
        self.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def start(self) -> None:
        sync_playwright = _import_playwright()
        self.profile_dir.mkdir(parents=True, exist_ok=True)
        self._ensure_writable_temp()
        self._pw = sync_playwright().start()
        try:
            self.ctx = self._launch()
            self.page = self.ctx.pages[0] if self.ctx.pages else self.ctx.new_page()
            self.page.set_default_timeout(20000)
            self.page.on("response", self._on_response)
            self._watched_pages.add(self.page)
            if self.page_router is not None:
                self.ctx.route("https://api.registry.igem.org/v1/blast/**", self.page_router)
            self._load_login_state()
        except Exception:
            # 启动失败也要收拾干净，否则下一次调用会报 asyncio loop 相关的怪错误
            self.close()
            raise

    def _login_state_path(self) -> Path:
        return self.profile_dir.parent / "login-state.json"

    def _load_login_state(self) -> None:
        """把随包携带的登录状态（login-state.json）塞进浏览器。"""
        path = self._login_state_path()
        if not path.exists():
            return
        try:
            state = json.loads(path.read_text(encoding="utf-8"))
            cookies = state.get("cookies") or []
            if cookies:
                self.ctx.add_cookies(cookies)
                self.log(f"已载入 login-state.json 里的 {len(cookies)} 个 cookie")
        except Exception as exc:  # noqa: BLE001
            self.log(f"载入登录状态失败：{exc}")

    def _save_login_state(self) -> None:
        """把当前登录状态存成 login-state.json，换机器/换浏览器也能直接用。"""
        try:
            state = self.ctx.storage_state()
            path = self._login_state_path()
            path.write_text(json.dumps(state, ensure_ascii=False), encoding="utf-8")
            self.log(f"登录状态已保存到 {path.name}")
        except Exception as exc:  # noqa: BLE001
            self.log(f"保存登录状态失败：{exc}")

    def _ensure_writable_temp(self) -> None:
        """Playwright 要在临时目录里建 artifacts；临时目录不可写时退到程序目录。"""
        import tempfile

        def writable(path: Optional[str]) -> bool:
            if not path:
                return False
            try:
                probe = Path(path) / f".igem_write_test_{os.getpid()}"
                probe.write_text("x", encoding="utf-8")
                probe.unlink()
                return True
            except Exception:  # noqa: BLE001
                return False

        extras = [os.environ.get("TEMP"), os.environ.get("TMP"),
                  os.environ.get("TMPDIR")]
        if any(writable(p) for p in extras):
            return

        local = self.profile_dir.parent / ".tmp"
        try:
            local.mkdir(parents=True, exist_ok=True)
            probe = local / f".igem_write_test_{os.getpid()}"
            probe.write_text("x", encoding="utf-8")
            probe.unlink()
            tempfile.tempdir = str(local)
            os.environ["TMP"] = os.environ["TEMP"] = os.environ["TMPDIR"] = str(local)
            self.log(f"临时目录不可写，改用 {local}")
        except Exception as exc:  # noqa: BLE001
            self.log(f"准备临时目录失败：{exc}")

    def _launch_once(self, channel):
        kwargs: Dict[str, Any] = {
            "user_data_dir": str(self.profile_dir),
            "headless": self.headless,
            "args": [
                "--no-first-run",
                "--no-default-browser-check",
                "--window-size=1380,900",
                "--window-position=60,60",
            ],
        }
        if self.headless:
            kwargs["viewport"] = {"width": 1366, "height": 900}
        else:
            # 让别人能看见窗口：不要用固定 viewport 去改窗口大小
            kwargs["no_viewport"] = True
        if channel:
            kwargs["channel"] = channel
        return self._pw.chromium.launch_persistent_context(**kwargs)

    def _close_stale_windows(self) -> bool:
        """关掉上次遗留的登录窗口：它占着配置目录会让新窗口起不来。"""
        if sys.platform != "win32":
            return False
        try:
            import ctypes
            from ctypes import wintypes

            user32 = ctypes.windll.user32
            kernel32 = ctypes.windll.kernel32
            pids = []
            proc = ctypes.WINFUNCTYPE(wintypes.BOOL, wintypes.HWND, wintypes.LPARAM)

            def callback(hwnd, _lparam):
                length = user32.GetWindowTextLengthW(hwnd)
                if length <= 0:
                    return True
                buf = ctypes.create_unicode_buffer(length + 1)
                user32.GetWindowTextW(hwnd, buf, length + 1)
                title = (buf.value or "").lower()
                if "chrome" in title and ("【igem 登录】" in title
                                          or "accounts.igem.org" in title):
                    pid = wintypes.DWORD()
                    user32.GetWindowThreadProcessId(hwnd, ctypes.byref(pid))
                    if pid.value:
                        pids.append(pid.value)
                return True

            user32.EnumWindows(proc(callback), 0)
            for pid in set(pids):
                handle = kernel32.OpenProcess(0x0001, False, pid)  # PROCESS_TERMINATE
                if handle:
                    kernel32.TerminateProcess(handle, 1)
                    kernel32.CloseHandle(handle)
                    self.log(f"已关闭残留的登录窗口（pid={pid}）")
            return bool(pids)
        except Exception as exc:  # noqa: BLE001
            self.log(f"关闭残留窗口失败：{exc}")
            return False

    def _launch(self):
        last_error: Optional[Exception] = None
        for attempt in (0, 1):
            for channel in ("chrome", "msedge", None):
                try:
                    return self._launch_once(channel)
                except Exception as exc:  # noqa: BLE001
                    last_error = exc
            if attempt == 0 and self._close_stale_windows():
                time.sleep(1.5)
                continue
            break
        raise RuntimeError(
            "无法启动浏览器（可能上一次的窗口还没关，或配置文件被占用）。\n"
            "如果任务栏里还有标题为【iGEM 登录】的 Chrome 窗口，先把它关掉。\n"
            f"配置文件目录：{self.profile_dir}\n"
            f"原始错误：{last_error}"
        )

    def close(self) -> None:
        try:
            if self.ctx is not None:
                self.ctx.close()
        except Exception:  # noqa: BLE001
            pass
        try:
            if self._pw is not None:
                self._pw.stop()
        except Exception:  # noqa: BLE001
            pass
        self.ctx = None
        self.page = None
        self._pw = None

    # -- 工具 -------------------------------------------------------------

    def log(self, msg: str) -> None:
        if self.verbose:
            print(f"[igem] {msg}", file=sys.stderr, flush=True)

    def _progress(self, msg: str) -> None:
        """把当前进度报给界面（如果有的话），同时写日志。"""
        if self.progress is not None:
            try:
                self.progress(msg)
            except Exception:  # noqa: BLE001
                pass
        self.log(msg)

    def _on_response(self, response) -> None:
        try:
            if "/v1/blast/jobs/" in response.url and response.request.method == "GET":
                self._job_responses.append(response)
        except Exception:  # noqa: BLE001
            pass

    def _me_status(self) -> Optional[int]:
        try:
            resp = self.ctx.request.get(ME_URL, fail_on_status_code=False, timeout=6000)
            return resp.status
        except Exception as exc:  # noqa: BLE001
            self.log(f"检查登录状态失败：{exc}")
            return None

    def _me_status_in_page(self) -> Optional[int]:
        """用页面自己的 fetch 再查一次（浏览器里的 cookie 才是最真实的）。"""
        if self.page is None or self.page.is_closed():
            return None
        try:
            status = self.page.evaluate(
                """async () => {
                    try {
                        const r = await fetch(
                            'https://api.registry.igem.org/v1/accounts/me',
                            { credentials: 'include' });
                        return r.status;
                    } catch (e) { return 0; }
                }"""
            )
            return int(status) if status else None
        except Exception as exc:  # noqa: BLE001
            self.log(f"页面内检查登录状态失败：{exc}")
            return None

    def _page_signed_in(self, page, wait_seconds: float = 2.5) -> bool:
        """看页面本身：有指向 /profile 的链接、又没有登录入口，就是已登录。

        比打接口可靠（接口有时对已登录会话返回 403/429），也比按文字匹配稳
        （改用 href 判断，不怕界面改文案）。
        """
        deadline = time.time() + wait_seconds
        while True:
            try:
                if page is None or page.is_closed():
                    return False
                url = page.url or ""
                if "accounts.igem.org" in url:
                    return False      # 还停在登录页
                if "registry.igem.org" not in url:
                    return False
                if page.locator('a[href="/profile"]').count() > 0:
                    return True
                has_sign_in = page.locator(
                    'a[href*="sign-in"], a[href*="signin"]').count() > 0
                if has_sign_in:
                    return False
                if page.locator('textarea[name="sequence"]').count() > 0:
                    # 已经在搜索页上、又没有登录入口，说明登录过了
                    return True
            except Exception as exc:  # noqa: BLE001
                self.log(f"页面登录状态检查失败：{exc}")
                return False
            if time.time() >= deadline:
                return False
            try:
                page.wait_for_timeout(200)
            except Exception:  # noqa: BLE001
                return False

    def _adopt_page(self, page) -> None:
        """把页面切换成某个标签页（比如登录是在新标签里完成的）。"""
        if page is self.page:
            return
        self.log("切换到最后打开的页面")
        self.page = page
        try:
            page.set_default_timeout(20000)
        except Exception:  # noqa: BLE001
            pass
        try:
            if page not in self._watched_pages:
                page.on("response", self._on_response)
                self._watched_pages.add(page)
        except Exception:  # noqa: BLE001
            pass

    def _ui_says_signed_in(self, wait_seconds: float = 2.5) -> bool:
        try:
            pages = [p for p in self.ctx.pages if not p.is_closed()]
        except Exception:  # noqa: BLE001
            pages = [self.page] if self.page is not None else []
        for page in pages:
            try:
                if "registry.igem.org" not in (page.url or ""):
                    continue
            except Exception:  # noqa: BLE001
                continue
            if self._page_signed_in(page, wait_seconds):
                self._adopt_page(page)
                return True
        return False

    def is_logged_in(self) -> bool:
        if self._ui_says_signed_in():
            return True
        try:
            on_registry = any(
                "registry.igem.org" in (p.url or "")
                for p in self.ctx.pages if not p.is_closed()
            )
        except Exception:  # noqa: BLE001
            on_registry = False
        if on_registry:
            return False
        return self._me_status() == 200

    def _save_debug(self, name: str) -> None:
        if self.debug_dir is None or self.page is None:
            return
        self.debug_dir.mkdir(parents=True, exist_ok=True)
        try:
            self.page.screenshot(path=str(self.debug_dir / f"{name}.png"), full_page=True)
            (self.debug_dir / f"{name}.html").write_text(
                self.page.content(), encoding="utf-8", errors="replace")
            self.log(f"调试文件已保存到 {self.debug_dir}")
        except Exception as exc:  # noqa: BLE001
            self.log(f"保存调试信息失败：{exc}")

    # -- 登录 -------------------------------------------------------------

    def _focus_browser_window(self) -> bool:
        """Windows 下把登录用的浏览器窗口提到最前，免得用户以为没反应。"""
        if sys.platform != "win32":
            return False
        try:
            import ctypes
            from ctypes import wintypes

            user32 = ctypes.windll.user32
            user32.SetWindowPos.argtypes = [
                wintypes.HWND, wintypes.HWND,
                ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                ctypes.c_uint,
            ]
            matches = []
            proc = ctypes.WINFUNCTYPE(wintypes.BOOL, wintypes.HWND, wintypes.LPARAM)

            def callback(hwnd, _lparam):
                if not user32.IsWindowVisible(hwnd):
                    return True
                length = user32.GetWindowTextLengthW(hwnd)
                if length <= 0:
                    return True
                buf = ctypes.create_unicode_buffer(length + 1)
                user32.GetWindowTextW(hwnd, buf, length + 1)
                title = (buf.value or "").lower()
                # 只认我们这个登录页的窗口，别把用户自己的浏览器窗口也抬上来
                if "chrome" in title and ("accounts.igem.org" in title
                                          or "registry.igem.org" in title
                                          or "【igem 登录】" in title):
                    matches.append(hwnd)
                return True

            user32.EnumWindows(proc(callback), 0)
            self.log(f"查找登录窗口：匹配到 {len(matches)} 个")
            if not matches:
                return False
            swp_nosize_nomove = 0x0001 | 0x0002

            class FLASHWINFO(ctypes.Structure):
                _fields_ = [
                    ("cbSize", wintypes.UINT),
                    ("hwnd", wintypes.HWND),
                    ("dwFlags", wintypes.DWORD),
                    ("uCount", wintypes.UINT),
                    ("dwTimeout", wintypes.DWORD),
                ]

            for hwnd in matches:
                user32.ShowWindow(hwnd, 9)          # SW_RESTORE
                # 直接置顶：登录是一次性操作，窗口关掉后置顶状态也随之消失
                user32.SetWindowPos(hwnd, wintypes.HWND(-1), 0, 0, 0, 0,
                                    swp_nosize_nomove)
                user32.SetForegroundWindow(hwnd)
                # 任务栏图标闪烁，提醒用户切过去
                info = FLASHWINFO(
                    ctypes.sizeof(FLASHWINFO), hwnd,
                    0x00000003 | 0x0000000C,  # FLASHW_ALL | FLASHW_TIMERNOFG
                    0, 0,
                )
                user32.FlashWindowEx(ctypes.byref(info))
            return True
        except Exception as exc:  # noqa: BLE001
            self.log(f"提升窗口失败：{exc}")
            return False

    def _ensure_window_visible(self) -> None:
        """把窗口恢复到屏幕内正常大小（浏览器可能记住上次最小化/跑出屏幕的位置）。"""
        if self.headless or self.page is None:
            return
        try:
            cdp = self.ctx.new_cdp_session(self.page)
            window = cdp.send("Browser.getWindowForTarget")
            cdp.send("Browser.setWindowBounds", {
                "windowId": window["windowId"],
                "bounds": {
                    "left": 60, "top": 60, "width": 1280, "height": 860,
                    "windowState": "normal",
                },
            })
            cdp.detach()
        except Exception as exc:  # noqa: BLE001
            self.log(f"调整窗口位置失败：{exc}")

    def ensure_login(self, interactive: bool) -> bool:
        """确认已登录；如果现场登录了，返回 True。"""
        if self.assume_logged_in:
            return False
        if self.is_logged_in():
            self.log("已登录。")
            return False

        if not interactive:
            raise LoginRequired(
                "尚未登录 iGEM Registry，BLAST 接口需要登录后才能使用。\n"
                "请先运行一次登录命令（会打开浏览器窗口，登录一次即可长期保存）：\n"
                f'    python "{Path(sys.argv[0]).name}" --login'
            )

        self._progress("等你登录 iGEM…")
        print("未检测到 iGEM 登录状态。", file=sys.stderr, flush=True)
        print("请在弹出的浏览器窗口里完成登录；", file=sys.stderr, flush=True)
        print("（如果没看到浏览器窗口，请点任务栏上的 Chrome 图标）", file=sys.stderr, flush=True)
        print(f"登录成功后程序会自动继续（最多等待 {int(self.login_timeout)} 秒）。",
              file=sys.stderr, flush=True)
        focused = False
        try:
            if "registry.igem.org" not in (self.page.url or ""):
                self.page.goto(SITE_URL, wait_until="domcontentloaded", timeout=60000)
            self.page.bring_to_front()
            # 直接跳到登录页，省得用户自己找 Sign In 按钮
            sign_in = ("https://accounts.igem.org/sign-in?redirect_to="
                       + quote(SITE_URL, safe=""))
            response = self.page.goto(sign_in, wait_until="domcontentloaded", timeout=60000)
            if response is not None and response.status >= 400:
                raise RuntimeError(f"sign-in page returned {response.status}")
            self.page.bring_to_front()
            self.page.wait_for_timeout(500)
            # 给这个窗口一个显眼的名字，方便用户在任务栏 / Alt+Tab 里认出来
            try:
                self.page.evaluate(
                    "document.title = '【iGEM 登录】请在这个窗口登录'")
            except Exception:  # noqa: BLE001
                pass
            self._ensure_window_visible()
            focused = self._focus_browser_window()
        except Exception:  # noqa: BLE001
            # 兜底：回到搜索页，点右上角的 Sign In
            try:
                self.page.goto(SITE_URL, wait_until="domcontentloaded", timeout=60000)
                self.page.get_by_role(
                    "link", name=re.compile(r"sign\s*in", re.I)
                ).first.click(timeout=10000)
                focused = self._focus_browser_window()
            except Exception:  # noqa: BLE001
                pass

        deadline = time.time() + self.login_timeout
        away_since = None
        while time.time() < deadline:
            if self.is_logged_in():
                print("登录成功。", file=sys.stderr, flush=True)
                self._save_login_state()
                self._progress("登录成功，继续查询…")
                return True
            # 兜底：登录页把用户送回 registry，就说明登录成功了；
            # 这里再稳几秒等页面渲染，避免导航栏还没出来就一直干等。
            try:
                url = self.page.url or ""
                if ("accounts.igem.org" not in url
                        and "registry.igem.org" in url):
                    if away_since is None:
                        away_since = time.time()
                    elif time.time() - away_since > 3:
                        print("已离开登录页，继续。", file=sys.stderr, flush=True)
                        self._save_login_state()
                        self._progress("登录完成，继续查询…")
                        return True
                else:
                    away_since = None
            except Exception:  # noqa: BLE001
                away_since = None
            if not focused:
                focused = self._focus_browser_window()
            try:
                if self.page.is_closed():
                    raise LoginRequired(
                        "浏览器窗口被关掉了，登录没有完成。请重新运行 login.bat。"
                    )
                self.page.wait_for_timeout(1000)
            except LoginRequired:
                raise
            except Exception:  # noqa: BLE001
                raise LoginRequired(
                    "浏览器窗口被关掉了，登录没有完成。请重新运行 login.bat。"
                )
        try:
            current_url = self.page.url
        except Exception:  # noqa: BLE001
            current_url = "?"
        raise LoginRequired(
            f"等待登录超时（{int(self.login_timeout)} 秒），一直没检测到登录。\n"
            "请确认你是在程序弹出的那个 Chrome 窗口里登录的——\n"
            "窗口标题是【iGEM 登录】；在你自己平时用的 Chrome 里登录不会共享给程序。\n"
            "登录成功的标志：那个窗口右上角变成你的账号，不再是 Sign In。\n"
            f"（当时页面：{current_url}）"
        )

    # -- 搜索 -------------------------------------------------------------

    def open_search_page(self) -> None:
        self._progress("打开 BLAST 页面…")
        self.page.goto(SITE_URL, wait_until="domcontentloaded", timeout=60000)
        self.page.wait_for_selector('textarea[name="sequence"]', timeout=30000)

    def submit(self, sequence: str) -> str:
        """填入序列并点击 Search，返回 jobId。"""
        self._job_responses = []

        if self.algorithm != DEFAULT_ALGORITHM:
            try:
                self.page.get_by_role("radio", name=self.algorithm, exact=False).click(timeout=5000)
            except Exception:  # noqa: BLE001
                # 页面上是 Radix 的 Tabs，不是原生 radio
                try:
                    self.page.get_by_role("tab", name=re.compile(self.algorithm, re.I)).click(timeout=5000)
                except Exception:  # noqa: BLE001
                    self.log(f"选择算法 {self.algorithm} 失败，使用默认算法。")

        self._progress(f"提交序列（{len(sequence)} bp）…")
        self.page.fill('textarea[name="sequence"]', sequence)
        self.page.get_by_role("button", name="Search", exact=True).click()

        try:
            self.page.wait_for_url(re.compile(r"/parts/blast/[0-9A-Za-z_.:-]+$"),
                                   timeout=60000)
        except Exception as exc:  # noqa: BLE001
            self._save_debug("submit_failed")
            if not self.assume_logged_in and not self.is_logged_in():
                raise LoginRequired(
                    "提交搜索时被拒绝（未登录）。请先运行 --login 完成登录。"
                ) from exc
            text = ""
            try:
                text = self.page.inner_text("body")[:400]
            except Exception:  # noqa: BLE001
                pass
            raise SearchFailed(f"点击 Search 后页面没有跳转到结果页。当前页面：{text!r}") from exc

        match = JOB_ID_RE.search(self.page.url)
        if not match:
            raise SearchFailed(f"无法从地址栏解析 jobId：{self.page.url}")
        job_id = match.group(1)
        self.log(f"已跳转到结果页，jobId = {job_id}")
        return job_id

    def wait_job(self, job_id: str) -> Dict[str, Any]:
        """等待页面自己轮询回来的 job 结果（不额外发请求）。"""
        deadline = time.time() + self.timeout_ms / 1000.0
        started = time.time()
        last_report = 0.0
        while time.time() < deadline:
            for response in reversed(self._job_responses):
                if job_id not in response.url:
                    continue
                try:
                    data = response.json()
                except Exception:  # noqa: BLE001
                    continue
                if not isinstance(data, dict):
                    continue
                if data.get("failedReason"):
                    raise SearchFailed(f"BLAST 任务失败：{data['failedReason']}")
                if data.get("finishedOn"):
                    if isinstance(data.get("returnvalue"), dict):
                        return data
                    raise SearchFailed(
                        "BLAST 任务已结束，但没有返回结果（可能队列里出错了）。"
                    )
            elapsed = time.time() - started
            if elapsed - last_report >= 2:
                last_report = elapsed
                self._progress(f"等 BLAST 出结果… 已等 {int(elapsed)} 秒")
            try:
                self.page.wait_for_timeout(1000)
            except Exception:  # noqa: BLE001
                time.sleep(1)
        raise SearchFailed(f"等待 BLAST 结果超时（{self.timeout_ms // 1000} 秒）。")

    def wait_results_rendered(self) -> None:
        """等到结果区域渲染完成（过滤按钮或空结果提示出现）。"""
        candidates = [
            self.page.get_by_role("button", name="Full match", exact=True),
            self.page.get_by_text("No matching sequences found."),
            self.page.get_by_text("No hits match the selected filters."),
        ]
        deadline = time.time() + 60
        while time.time() < deadline:
            for loc in candidates:
                try:
                    if loc.count() > 0 and loc.first.is_visible():
                        return
                except Exception:  # noqa: BLE001
                    pass
            self.page.wait_for_timeout(500)
        self._save_debug("results_not_rendered")
        raise SearchFailed("结果页没有渲染出结果列表。")

    def apply_full_match_filter(self) -> bool:
        """点击 "Full match" 过滤按钮。返回是否成功点中。"""
        button = self.page.get_by_role("button", name="Full match", exact=True)
        if button.count() == 0:
            self.log("页面没有 Full match 按钮（可能 0 条命中）。")
            return False
        try:
            button.first.click(timeout=10000)
        except Exception as exc:  # noqa: BLE001
            self.log(f"点击 Full match 失败：{exc}")
            return False
        self.page.wait_for_timeout(600)
        try:
            pressed = button.first.get_attribute("aria-pressed")
            self.log(f"Full match 过滤状态：aria-pressed={pressed}")
        except Exception:  # noqa: BLE001
            pass
        return True

    def read_first_ui_result(self):
        """读取 Full match 过滤后第一条结果的 (砖号, 卡片文本)。"""
        cards = self.page.locator("div.rounded-xl.border.bg-white.p-5")
        try:
            if cards.count() > 0:
                card = cards.first
                link = card.locator('a[href^="/parts/"]').first
                if link.count() > 0:
                    name = (link.inner_text() or "").strip()
                    text = card.inner_text() or ""
                    if name:
                        return name, text
        except Exception as exc:  # noqa: BLE001
            self.log(f"读取结果卡片失败：{exc}")

        # 兜底：扫描结果区的部件链接
        links = self.page.locator('a[href^="/parts/"]')
        for index in range(links.count()):
            link = links.nth(index)
            try:
                href = (link.get_attribute("href") or "").rstrip("/")
                if href in ("/parts/blast", "/parts"):
                    continue
                text = (link.inner_text() or "").strip()
                if text:
                    return text, ""
            except Exception:  # noqa: BLE001
                continue
        return None, ""

    def page_says_no_filtered_hits(self) -> bool:
        try:
            return self.page.get_by_text("No hits match the selected filters.").count() > 0
        except Exception:  # noqa: BLE001
            return False

    # -- 主流程 -----------------------------------------------------------

    def run(self, sequence: str, *, interactive_login: bool = False) -> Result:
        self._progress("启动浏览器…")
        self.open_search_page()
        did_login = self.ensure_login(interactive=interactive_login)
        if did_login:
            # 登录后一般会自动跳回搜索页；只有在表单不在时才重新加载
            try:
                if self.page.locator('textarea[name="sequence"]').count() == 0:
                    self.open_search_page()
            except Exception:  # noqa: BLE001
                self.open_search_page()

        job_id = self.submit(sequence)
        job = self.wait_job(job_id)

        data = job.get("data") or {}
        job_sequence = clean_sequence(str(data.get("sequence") or ""))
        query_length = len(job_sequence) or len(sequence)

        return_value = job.get("returnvalue") or {}
        if return_value.get("success") is False:
            raise SearchFailed(str(return_value.get("error") or "BLAST 返回失败。"))

        raw_hits = return_value.get("hits") or []
        hits = [Hit.from_api(item) for item in raw_hits if isinstance(item, dict)]

        self.wait_results_rendered()
        self.apply_full_match_filter()
        ui_name, _ui_text = self.read_first_ui_result()
        if self.page_says_no_filtered_hits():
            ui_name = None

        verdict = judge(hits, query_length)
        first: Optional[Hit] = verdict["first"]

        # 以页面显示的第一条为准做交叉校验
        if first is not None and ui_name and ui_name != first.part_name:
            self.log(f"页面首条为 {ui_name}，接口首条为 {first.part_name}，以接口为准。")

        self._save_debug("result")

        return Result(
            status=verdict["status"],
            part=(first.part_name if first is not None else None),
            query_length=query_length,
            identity=(first.identity if first is not None else None),
            subject_length=(first.subject_length if first is not None else None),
            job_id=job_id,
            url=self.page.url,
            message=verdict["message"],
            first_hit=first,
            hits_total=len(hits),
            hits_full_match=verdict["full_count"],
            source="api",
        )


# --------------------------------------------------------------------------
# 对外入口
# --------------------------------------------------------------------------

def lookup(
    sequence: str,
    *,
    headless: bool = True,
    profile_dir: Optional[Path] = None,
    algorithm: str = DEFAULT_ALGORITHM,
    timeout: float = 180.0,
    login_timeout: float = 300.0,
    verbose: bool = False,
    debug_dir: Optional[Path] = None,
    interactive_login: bool = False,
    progress: Optional[Callable[[str], None]] = None,
    _router: Optional[Callable[[Any], None]] = None,
    _assume_logged_in: bool = False,
) -> Result:
    """输入 DNA 序列，返回查询结果对象。"""
    seq = validate_sequence(clean_sequence(sequence))
    with Session(
        headless=headless,
        profile_dir=profile_dir,
        algorithm=algorithm,
        timeout=timeout,
        login_timeout=login_timeout,
        verbose=verbose,
        debug_dir=debug_dir,
        page_router=_router,
        assume_logged_in=_assume_logged_in,
        progress=progress,
    ) as session:
        return session.run(seq, interactive_login=interactive_login)


def login(profile_dir: Optional[Path] = None, timeout: float = 300.0,
          verbose: bool = False) -> None:
    """打开浏览器窗口，等用户完成登录后保存登录态。"""
    with Session(headless=False, profile_dir=profile_dir, verbose=verbose,
                 login_timeout=timeout) as session:
        session.page.goto(SITE_URL, wait_until="domcontentloaded", timeout=60000)
        if session.is_logged_in():
            print("已经处于登录状态，无需重复登录。", file=sys.stderr)
            return
        session.ensure_login(interactive=True)


# --------------------------------------------------------------------------
# 自检（不需要账号：用假的 BLAST 接口驱动真实页面）
# --------------------------------------------------------------------------

MOCK_SEQUENCE = (
    "ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCA"
    "GTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACA"
    "CTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCGAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGC"
    "CCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATAC"
    "CCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCA"
    "CACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTC"
    "AACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGT"
    "CCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAA"
    "CTATACAAATAA"
)


def _mock_hit(part: str, identity: float, q_start: int, q_end: int,
              s_start: int, s_end: int, s_len: int,
              mismatches: int = 0, gaps: int = 0) -> Dict[str, Any]:
    return {
        "partName": part, "identity": identity,
        "queryStart": q_start, "queryEnd": q_end,
        "subjectStart": s_start, "subjectEnd": s_end, "subjectLength": s_len,
        "alignmentLength": q_end - q_start + 1,
        "mismatches": mismatches, "gaps": gaps, "gapOpenings": 0,
        "bitScore": 900.0, "evalue": 0.0,
    }


def _mock_router_factory(job: Dict[str, Any]):
    def router(route) -> None:
        url = route.request.url
        if url.endswith("/blast/search"):
            route.fulfill(status=200, content_type="application/json",
                          body=json.dumps({"jobID": job["id"]}))
        elif "/blast/jobs/" in url:
            route.fulfill(status=200, content_type="application/json",
                          body=json.dumps(job))
        else:
            route.continue_()

    return router


def _mock_job(seq: str, hits: List[Dict[str, Any]]) -> Dict[str, Any]:
    return {
        "id": "mockjob123",
        "data": {"sequence": seq, "algorithm": "megablast"},
        "finishedOn": int(time.time() * 1000),
        "failedReason": None,
        "returnvalue": {"success": True, "searchTime": 12.5, "hits": hits},
    }


def self_test(verbose: bool = True) -> int:
    """用假接口跑一遍完整 UI 流程，验证程序本身是否正常。"""
    seq = MOCK_SEQUENCE
    n = len(seq)
    failures = 0

    cases = [
        (
            "完全一致的砖（排在 partial 命中之后）",
            [  # 96% 的局部命中会被 Full match 过滤掉
                _mock_hit("BBa_J04450", 96.0, 1, n, 20, n + 19, n + 19),
                _mock_hit("BBa_E0040", 100.0, 1, n, 1, n, n),
                _mock_hit("BBa_K1234567", 99.7, 1, n, 1, n - 5, n - 5, mismatches=1),
            ],
            "exact", "BBa_E0040",
        ),
        (
            "首条 Full match 只有 99.7%",
            [
                _mock_hit("BBa_K1234567", 99.7, 1, n, 1, n - 5, n - 5, mismatches=1),
            ],
            "not_exact", "BBa_K1234567",
        ),
        (
            "完全没有命中",
            [],
            "no_hits", None,
        ),
    ]

    for title, hits, want_status, want_part in cases:
        job = _mock_job(seq, hits)
        print(f"· 用例：{title}", file=sys.stderr, flush=True)
        try:
            result = lookup(seq, headless=True, verbose=verbose,
                            _router=_mock_router_factory(job),
                            _assume_logged_in=True)
        except Exception as exc:  # noqa: BLE001
            failures += 1
            print(f"  ✗ 异常：{exc}", file=sys.stderr)
            if verbose:
                traceback.print_exc()
            continue
        ok = result.status == want_status and result.part == want_part
        if ok:
            print(f"  ✓ status={result.status} part={result.part}", file=sys.stderr)
        else:
            failures += 1
            print(f"  ✗ 期望 status={want_status} part={want_part}，"
                  f"实际 status={result.status} part={result.part}", file=sys.stderr)
    return 1 if failures else 0


# --------------------------------------------------------------------------
# 图形界面
# --------------------------------------------------------------------------

def set_clipboard(root, text: str) -> bool:
    """把文本写进剪贴板。

    Windows 下直接用 API 写，这样程序关掉以后内容还在（Tk 自带的剪贴板
    在窗口销毁后可能就没了）。
    """
    if sys.platform == "win32":
        try:
            import ctypes
            from ctypes import wintypes

            user32 = ctypes.windll.user32
            kernel32 = ctypes.windll.kernel32
            CF_UNICODETEXT = 13
            GMEM_MOVEABLE = 0x0002

            kernel32.GlobalAlloc.argtypes = [wintypes.UINT, ctypes.c_size_t]
            kernel32.GlobalAlloc.restype = wintypes.HANDLE
            kernel32.GlobalLock.argtypes = [wintypes.HANDLE]
            kernel32.GlobalLock.restype = wintypes.LPVOID
            kernel32.GlobalUnlock.argtypes = [wintypes.HANDLE]
            kernel32.GlobalFree.argtypes = [wintypes.HANDLE]
            user32.OpenClipboard.argtypes = [wintypes.HWND]
            user32.SetClipboardData.argtypes = [wintypes.UINT, wintypes.HANDLE]
            user32.SetClipboardData.restype = wintypes.HANDLE

            if not user32.OpenClipboard(None):
                return False
            try:
                user32.EmptyClipboard()
                data = text.encode("utf-16-le") + b"\x00\x00"
                handle = kernel32.GlobalAlloc(GMEM_MOVEABLE, len(data))
                if not handle:
                    return False
                ptr = kernel32.GlobalLock(handle)
                if not ptr:
                    kernel32.GlobalFree(handle)
                    return False
                ctypes.memmove(ptr, data, len(data))
                kernel32.GlobalUnlock(handle)
                if not user32.SetClipboardData(CF_UNICODETEXT, handle):
                    kernel32.GlobalFree(handle)
                    return False
                return True
            finally:
                user32.CloseClipboard()
        except Exception:  # noqa: BLE001
            pass
    try:
        root.clipboard_clear()
        root.clipboard_append(text)
        root.update()
        return True
    except Exception:  # noqa: BLE001
        return False


def run_gui(profile_dir: Optional[Path] = None, headless: bool = True,
            selftest: bool = False) -> None:
    import tkinter as tk
    import queue
    from tkinter import messagebox, ttk

    root = tk.Tk()
    root.title(f"iGEM 砖号查询器 v{APP_VERSION}")
    root.geometry("800x640")
    root.minsize(620, 480)

    style = ttk.Style()
    try:
        style.theme_use("vista")
    except tk.TclError:
        pass

    header = ttk.Frame(root, padding=(16, 14, 16, 6))
    header.pack(fill="x")
    ttk.Label(header, text="iGEM 砖号查询",
              font=("Segoe UI", 16, "bold")).pack(anchor="w")
    ttk.Label(header, foreground="#555", wraplength=740,
              text="粘贴 DNA 序列 → 自动在 registry.igem.org/parts/blast 上搜索 → "
                   "点 Full match → 看第一条").pack(anchor="w", pady=(2, 0))

    body = ttk.Frame(root, padding=(16, 6, 16, 6))
    body.pack(fill="both", expand=True)

    ttk.Label(body, text="DNA 序列（可直接粘贴 FASTA）").pack(anchor="w")
    text_frame = ttk.Frame(body)
    text_frame.pack(fill="both", expand=True, pady=(4, 8))
    text = tk.Text(text_frame, height=10, wrap="char", font=("Consolas", 10),
                   relief="solid", borderwidth=1, undo=True)
    scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=text.yview)
    text.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side="right", fill="y")
    text.pack(side="left", fill="both", expand=True)

    buttons = ttk.Frame(body)
    buttons.pack(fill="x", pady=(0, 8))

    search_btn = ttk.Button(buttons, text="查询砖号")
    search_btn.pack(side="left")
    login_btn = ttk.Button(buttons, text="登录 iGEM")
    login_btn.pack(side="left", padx=(8, 0))
    clear_btn = ttk.Button(buttons, text="清空")
    clear_btn.pack(side="left", padx=(8, 0))

    result_row = ttk.Frame(body)
    result_row.pack(fill="x", pady=(4, 2))
    result_var = tk.StringVar(value="等待粘贴序列…")
    result_label = tk.Label(result_row, textvariable=result_var,
                            font=("Consolas", 20, "bold"),
                            fg="#1a7f37", anchor="w", justify="left")
    result_label.pack(side="left")
    copy_btn = ttk.Button(result_row, text="复制砖号", width=11, state="disabled")
    copy_btn.pack(side="left", padx=(14, 0))

    detail_var = tk.StringVar(
        value="第 1 步：把 DNA 序列粘到上面的框里　第 2 步：点「查询砖号」")
    detail = ttk.Label(body, textvariable=detail_var, foreground="#555",
                       wraplength=700, justify="left")
    detail.pack(fill="x")

    status_var = tk.StringVar()
    shared_state = (Path(profile_dir) if profile_dir else DEFAULT_PROFILE).parent / "login-state.json"
    status_var.set(f"就绪 · v{APP_VERSION} · "
                   + ("共享登录已载入" if shared_state.exists() else "未配置共享登录"))
    status = ttk.Label(root, textvariable=status_var, foreground="#666",
                       padding=(16, 4, 16, 8))
    status.pack(fill="x", side="bottom")

    state = {"busy": False, "last_seq": "", "part": ""}
    events: "queue.Queue" = queue.Queue()

    def set_busy(busy: bool) -> None:
        state["busy"] = busy
        st = "disabled" if busy else "normal"
        search_btn.configure(state=st)
        login_btn.configure(state=st)
        clear_btn.configure(state=st)

    def on_copy() -> None:
        part = state.get("part") or ""
        if not part:
            return
        if set_clipboard(root, part):
            copy_btn.configure(text="已复制 ✓")
            status_var.set(f"已复制 {part} —— 可以直接 Ctrl+V 粘贴")
            root.after(1500, lambda: copy_btn.configure(text="复制砖号"))
        else:
            status_var.set("复制失败，请手动选中砖号复制")

    def finish(result) -> None:
        set_busy(False)
        if result.status == "exact":
            state["part"] = result.part or ""
            result_var.set(result.part or "")
            result_label.configure(fg="#1a7f37")
            copy_btn.configure(state="normal", text="复制砖号")
            detail_var.set(f"{result.message}\n"
                           f"job: {result.job_id} · 命中 {result.hits_total} 条，"
                           f"其中 Full match {result.hits_full_match} 条")
        elif result.status == "not_exact":
            state["part"] = ""
            copy_btn.configure(state="disabled", text="复制砖号")
            result_var.set("不是完全一致的砖")
            result_label.configure(fg="#b35900")
            detail_var.set(result.message)
        else:
            state["part"] = ""
            copy_btn.configure(state="disabled", text="复制砖号")
            result_var.set("没有找到")
            result_label.configure(fg="#b35900")
            detail_var.set(result.message)
        status_var.set("完成")

    def fail(err: Exception) -> None:
        set_busy(False)
        if isinstance(err, LoginRequired):
            result_var.set("需要先登录 iGEM")
            result_label.configure(fg="#b35900")
            detail_var.set("Registry 的 BLAST 必须登录才能用。\n"
                           "点「是」会打开一个 Chrome 窗口（标题带【iGEM 登录】），\n"
                           "在里面登录一次；登录成功后会自动接着查询，不用重新点。")
            status_var.set("未登录")
            if messagebox.askyesno(
                "需要登录",
                "还没有登录 iGEM。\n\n现在打开浏览器登录吗？\n"
                "（登录一次即可长期保存）",
            ):
                start_query(state.get("last_seq", ""), interactive=True)
        else:
            result_var.set("出错了")
            result_label.configure(fg="#b02a37")
            detail_var.set(str(err))
            status_var.set("出错")

    def worker(seq: str, interactive: bool = False) -> None:
        try:
            res = lookup(
                seq,
                headless=(not interactive) and headless,
                profile_dir=profile_dir,
                interactive_login=interactive,
                progress=lambda msg: events.put(("status", msg)),
            )
            events.put(("ok", res))
        except Exception as exc:  # noqa: BLE001
            events.put(("err", exc))

    def start_query(clean: str, interactive: bool = False) -> None:
        if not clean:
            return
        state["last_seq"] = clean
        state["part"] = ""
        copy_btn.configure(state="disabled", text="复制砖号")
        set_busy(True)
        if interactive:
            result_var.set("等你登录…")
            result_label.configure(fg="#b35900")
            detail_var.set("登录窗口马上弹出来（标题【iGEM 登录】）。\n"
                           "登录成功后会自动继续查询，不用再点。\n"
                           "没看到窗口就去任务栏点 Chrome。")
            status_var.set("正在打开登录窗口…")
        else:
            result_var.set("查询中…")
            result_label.configure(fg="#333")
            detail_var.set(f"序列长度 {len(clean)} bp")
            status_var.set("正在启动浏览器…")
        threading.Thread(target=worker, args=(clean, interactive), daemon=True).start()

    def on_search() -> None:
        if state["busy"]:
            return
        seq = text.get("1.0", "end")
        try:
            clean = validate_sequence(clean_sequence(seq))
        except ValueError as exc:
            result_var.set("还没法查询")
            result_label.configure(fg="#b35900")
            detail_var.set(f"{exc}\n请先把序列粘到上面的框里（Ctrl+V），再点「查询砖号」。")
            status_var.set("等待序列")
            text.focus_set()
            return
        start_query(clean, interactive=False)

    def on_text_change(_event=None) -> None:
        """输入框内容变化时给点反馈，别让人对着「等待输入」发呆。"""
        try:
            if not text.edit_modified():
                return
            text.edit_modified(False)
        except Exception:  # noqa: BLE001
            return
        if state["busy"]:
            return
        seq = clean_sequence(text.get("1.0", "end"))
        if seq:
            status_var.set(f"已输入 {len(seq)} bp —— 点「查询砖号」开始（Ctrl+Enter 也行）")
            if result_var.get().startswith(("等待", "还没法", "请先")):
                result_var.set("准备就绪")
                result_label.configure(fg="#333")
                detail_var.set(f"序列长度 {len(seq)} bp　点「查询砖号」开始查询")
        else:
            status_var.set("就绪")
            result_var.set("等待粘贴序列…")
            result_label.configure(fg="#1a7f37")
            detail_var.set("第 1 步：把 DNA 序列粘到上面的框里　第 2 步：点「查询砖号」")

    def do_login() -> None:
        if state["busy"]:
            return
        set_busy(True)
        result_var.set("等你登录…")
        result_label.configure(fg="#b35900")
        detail_var.set("马上会弹出浏览器并打开 iGEM 登录页；\n"
                       "登录成功后浏览器会自动关闭，这里会自己继续。\n"
                       "如果没看到浏览器窗口，去任务栏点一下 Chrome 图标。")
        status_var.set("正在打开登录页面…")
        started = time.time()

        def tick() -> None:
            if not state["busy"]:
                return
            secs = int(time.time() - started)
            status_var.set(f"等待登录… 已等 {secs} 秒"
                           "（要在标题为【iGEM 登录】的 Chrome 窗口里登录）")
            root.after(2000, tick)

        root.after(2000, tick)

        def run() -> None:
            try:
                login(profile_dir=profile_dir)
                events.put(("login_ok", None))
            except Exception as exc:  # noqa: BLE001
                events.put(("err", exc))

        threading.Thread(target=run, daemon=True).start()

    # 后台线程只往队列里投递结果，UI 更新统一回到主线程做（Tk 不是线程安全的）
    def drain() -> None:
        try:
            while True:
                kind, payload = events.get_nowait()
                if kind == "ok":
                    finish(payload)
                elif kind == "err":
                    fail(payload)
                elif kind == "status":
                    status_var.set(payload)
                elif kind == "login_ok":
                    set_busy(False)
                    status_var.set("登录完成")
                    result_var.set("已登录")
                    result_label.configure(fg="#1a7f37")
                    detail_var.set("现在可以粘贴序列查询了。")
        except queue.Empty:
            pass
        root.after(120, drain)

    def on_clear() -> None:
        if state["busy"]:
            return
        text.delete("1.0", "end")
        state["part"] = ""
        copy_btn.configure(state="disabled", text="复制砖号")
        result_var.set("等待粘贴序列…")
        result_label.configure(fg="#1a7f37")
        detail_var.set("第 1 步：把 DNA 序列粘到上面的框里　第 2 步：点「查询砖号」")
        status_var.set("就绪")
        text.focus_set()

    search_btn.configure(command=on_search)
    login_btn.configure(command=do_login)
    clear_btn.configure(command=on_clear)
    copy_btn.configure(command=on_copy)
    result_label.bind("<Button-1>", lambda _e: on_copy())
    text.bind("<<Modified>>", on_text_change)
    root.bind("<Control-Return>", lambda _e: on_search())
    text.focus_set()
    root.after(120, drain)

    if selftest:
        # 隐藏的自检：塞一个假结果，点一下「复制砖号」，看剪贴板对不对
        finish(Result(status="exact", part="BBa_SELFTEST1", query_length=35,
                      identity=100.0, subject_length=35, job_id="selftest",
                      message="自检"))
        root.update()
        copy_btn.invoke()
        root.update()
        try:
            got = root.clipboard_get()
        except Exception:  # noqa: BLE001
            got = ""
        print("GUI_SELFTEST_CLIPBOARD=" + got, flush=True)
        root.destroy()
        return

    root.mainloop()


# --------------------------------------------------------------------------
# 命令行
# --------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="igem_brick_finder",
        description="输入 DNA 序列，在 iGEM Registry 的 BLAST 页面上搜索，"
                    "输出完全一致（100% 且等长）的砖号。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python igem_brick_finder.py --login                 # 首次：登录一次\n"
            '  python igem_brick_finder.py "ATGCGT..."             # 查询\n'
            "  python igem_brick_finder.py -f seq.fasta --json     # 从文件读，输出 JSON\n"
            "  python igem_brick_finder.py --gui                   # 打开小窗口\n"
            "  python igem_brick_finder.py --self-test             # 不需要账号的自检\n"
        ),
    )
    parser.add_argument("sequence", nargs="?", help="DNA 序列；写 - 表示从标准输入读取")
    parser.add_argument("-f", "--file", help="从文件读取序列（FASTA 或纯序列）")
    parser.add_argument("--login", action="store_true", help="只打开浏览器完成登录")
    parser.add_argument("--gui", action="store_true", help="打开图形界面")
    parser.add_argument("--gui-selftest", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--self-test", action="store_true",
                        help="用假数据跑通整个流程做自检（无需登录）")
    parser.add_argument("--headless", action="store_true",
                        help="不显示浏览器窗口（需要之前已登录过）")
    parser.add_argument("--profile", help=f"登录态保存目录（默认 {DEFAULT_PROFILE}）")
    parser.add_argument("--algorithm", choices=ALGORITHMS, default=DEFAULT_ALGORITHM,
                        help="BLAST 算法，默认 megablast")
    parser.add_argument("--timeout", type=float, default=180.0,
                        help="等待 BLAST 结果的秒数，默认 180")
    parser.add_argument("--login-timeout", type=float, default=300.0,
                        help="等待人工登录的秒数，默认 300")
    parser.add_argument("--json", action="store_true", help="输出 JSON")
    parser.add_argument("--out", help="把 JSON 结果写到指定文件（打包版没有控制台时用）")
    parser.add_argument("-v", "--verbose", action="store_true", help="打印过程信息")
    parser.add_argument("--debug", action="store_true",
                        help="保存截图/HTML 到 ./debug 便于排查")
    parser.add_argument("-q", "--quiet", action="store_true", help="只输出结果本身")
    return parser


def read_sequence(args) -> str:
    if args.file:
        return Path(args.file).read_text(encoding="utf-8", errors="replace")
    if args.sequence == "-":
        return sys.stdin.read()
    if args.sequence:
        return args.sequence
    if not sys.stdin.isatty():
        return sys.stdin.read()
    return input("请输入 DNA 序列: ")


def main(argv: Optional[List[str]] = None) -> int:
    # 打包成无控制台的 exe 时，stdout/stderr 是 None，先兜住免得 print 崩掉
    if sys.stdout is None:
        sys.stdout = open(os.devnull, "w", encoding="utf-8")
    if sys.stderr is None:
        sys.stderr = open(os.devnull, "w", encoding="utf-8")
    # Windows 控制台默认 GBK，统一成 UTF-8 以免中文/特殊字符崩掉
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass

    parser = build_parser()
    args = parser.parse_args(argv)
    profile = Path(args.profile) if args.profile else None

    if args.self_test:
        return self_test(verbose=args.verbose)

    if args.gui:
        run_gui(profile_dir=profile, headless=True)
        return 0

    if args.gui_selftest:
        run_gui(selftest=True)
        return 0

    if args.login:
        try:
            login(profile_dir=profile, timeout=args.login_timeout, verbose=args.verbose)
            return 0
        except Exception as exc:  # noqa: BLE001
            print(f"登录失败：{exc}", file=sys.stderr)
            return 1

    try:
        raw = read_sequence(args)
    except Exception as exc:  # noqa: BLE001
        print(f"读取序列失败：{exc}", file=sys.stderr)
        return 1

    try:
        result = lookup(
            raw,
            headless=args.headless,
            profile_dir=profile,
            algorithm=args.algorithm,
            timeout=args.timeout,
            login_timeout=args.login_timeout,
            verbose=args.verbose,
            debug_dir=(HERE / "debug") if args.debug else None,
            interactive_login=not args.headless,
        )
    except Exception as exc:  # noqa: BLE001
        if isinstance(exc, LoginRequired):
            message = str(exc)
        else:
            message = f"查询失败：{exc}"
        print(message, file=sys.stderr)
        if args.verbose or args.debug:
            traceback.print_exc()
        if args.out:
            try:
                Path(args.out).write_text(
                    json.dumps({"status": "error", "message": message},
                               ensure_ascii=False, indent=2),
                    encoding="utf-8")
            except Exception:  # noqa: BLE001
                pass
        return 1

    if args.out:
        try:
            Path(args.out).write_text(
                json.dumps(result.to_dict(), ensure_ascii=False, indent=2),
                encoding="utf-8")
        except Exception as exc:  # noqa: BLE001
            print(f"写结果文件失败：{exc}", file=sys.stderr)
            return 1
    elif args.json:
        print(json.dumps(result.to_dict(), ensure_ascii=False, indent=2))
    elif result.status == "exact":
        print(result.part)
        if not args.quiet:
            print(result.message, file=sys.stderr)
    else:
        print(result.message, file=sys.stderr)

    return 0 if result.status == "exact" else 2


if __name__ == "__main__":
    raise SystemExit(main())
