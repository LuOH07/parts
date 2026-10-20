#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
iGEM 砖号查询 · 网页版后端

一个文件同时干两件事：
  1. 把 index.html 发出去（所以本地运行只要打开 http://127.0.0.1:8765 就行）
  2. /api/convert：收 .dna + 用户自己的 iGEM 账号密码 → 登录 → 查砖号 → 把结果流式吐回来

为什么必须有这个后端？
  iGEM 的接口只允许 registry.igem.org 这个来源访问（跨域请求不返回
  Access-Control-Allow-Origin），浏览器里的网页直接调会被拦掉，
  所以登录和查询这两步只能放在服务端做。

安全：
  · 密码只用于这一次请求里的登录，不写日志、不落盘、不保存；
  · 每个请求自己开一个临时浏览器配置目录，用完就删；
  · 登录态不会在请求之间保留，所以也不存在“串号”。

运行：
  python app.py            # 默认 127.0.0.1:8765
  PORT=8000 python app.py  # 换端口（云平台会自己给 PORT）
"""

import base64
import json
import os
import shutil
import sys
import tempfile
import threading
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "app"))

import dna_to_bricks as engine          # noqa: E402
import igem_brick_finder as igem        # noqa: E402

PORT = int(os.environ.get("PORT", "8765"))
HOST = os.environ.get("HOST", "0.0.0.0")
ALLOW_ORIGIN = os.environ.get("ALLOW_ORIGIN", "*")      # 部署到别处时按需收紧
VERBOSE = os.environ.get("VERBOSE", "") == "1"

SIGN_IN = ("https://accounts.igem.org/sign-in?redirect_to="
           "https%3A%2F%2Fregistry.igem.org%2Fparts%2Fblast")

HEADLESS = os.environ.get("HEADLESS", "1") != "0"
USER_AGENT = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
              "(KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36")


def _install_stealth_launch() -> None:
    """让服务端的浏览器别一眼就是机器人：正常的 UA、语言、时区、窗口大小。

    云主机的 headless Chrome 很容易被挡；这里至少把能改的都改成普通 Chrome 的样子。
    """
    def _launch_once(self, channel):
        kwargs = {
            "user_data_dir": str(self.profile_dir),
            "headless": self.headless,
            "args": [
                "--no-first-run",
                "--no-default-browser-check",
                "--window-size=1366,900",
                "--window-position=40,40",
                "--disable-blink-features=AutomationControlled",
                "--lang=en-US",
            ],
            "user_agent": USER_AGENT,
            "locale": "en-US",
            "timezone_id": "Asia/Shanghai",
            "viewport": {"width": 1366, "height": 900},
        }
        if channel:
            kwargs["channel"] = channel
        return self._pw.chromium.launch_persistent_context(**kwargs)

    igem.Session._launch_once = _launch_once       # type: ignore[assignment]


_install_stealth_launch()

# 同一个进程里同时跑多个 Playwright 会打架，串起来执行
_RUN_LOCK = threading.Lock()


class CredError(RuntimeError):
    """账号密码不对 / 登录没成功。"""


# --------------------------------------------------------------------------
# 登录
# --------------------------------------------------------------------------

_NOT_AN_ERROR = {
    "username or email", "password", "sign in", "sign-in", "sign up!",
    "forgot password?", "signin",
}

# iGEM 英文报错 → 给人看的中文提示
_ALERT_HINTS = {
    "invalid credentials":
        "账号或密码不对。注意：如果这个 iGEM 账号平时是用 Google / 学校邮箱"
        "一键登录的，需要先在 iGEM 网站上单独设一个密码才能用这里。",
    "too many":
        "失败次数太多被暂时限制了，过一会儿再试。",
}


def _visible_alert(page) -> str:
    """把登录页上真正显示出来的错误提示抠出来。

    iGEM 的登录页把报错放在 [class*="danger"] 里（例如 Invalid credentials），
    而 [class*="invalid"] 会匹配到输入框外面的包装元素，里面是字段名，
    所以这里优先 danger，并且把字段名/按钮文字排除掉。
    """
    for sel in ('[class*="danger"]', '[role="alert"]',
                '[class*="error"]', '[class*="Error"]'):
        try:
            loc = page.locator(sel)
            n = loc.count()
        except Exception:  # noqa: BLE001
            continue
        for i in range(min(n, 6)):
            try:
                text = (loc.nth(i).inner_text() or "").strip()
            except Exception:  # noqa: BLE001
                continue
            text = " ".join(text.split())
            if not (2 < len(text) < 300):
                continue
            if text.lower() in _NOT_AN_ERROR:
                continue
            return text
    return ""


def _page_hint(page) -> str:
    """猜一下页面现在停在哪一步：验证码 / 两步验证 / 又回到登录框。"""
    try:
        if page.locator('iframe[src*="challenge"], iframe[src*="captcha"]').count():
            return "页面出现了人机验证（CAPTCHA），需要人工过一次"
        txt = (page.inner_text("body") or "").lower()
        if "verify you are human" in txt or "captcha" in txt:
            return "页面出现了人机验证（CAPTCHA）"
        for sel in ('input[name*="code" i]', 'input[autocomplete="one-time-code"]',
                    'input[name*="otp" i]', 'input[name*="token" i]'):
            if page.locator(sel).count():
                return "这个账号需要两步验证（要输入验证码），程序填不了"
        if page.locator("#password").count():
            return "还停在登录框（账号或密码被拒绝）"
    except Exception:  # noqa: BLE001
        pass
    return ""


def login_with_state(state_bytes: bytes, log=None) -> tuple:
    """用桌面版导出的 login-state.json（已登录的 cookie）建立会话，不走密码。"""
    say = log or (lambda _m: None)
    try:
        state = json.loads(state_bytes.decode("utf-8"))
    except Exception as exc:  # noqa: BLE001
        raise CredError("这个 json 解析不了，确认是桌面版 app 目录里的 "
                        "login-state.json。") from exc
    cookies = state.get("cookies") or []
    if not cookies:
        raise CredError("这个文件里没有 cookie，不是有效的登录状态。")

    profile = Path(tempfile.mkdtemp(prefix="igem-web-"))
    session = igem.Session(headless=HEADLESS, profile_dir=profile, verbose=False,
                           timeout=240.0)
    session.start()
    try:
        say(f"载入登录状态：{len(cookies)} 个 cookie")
        session.ctx.add_cookies(cookies)
        session.open_search_page()
        if not session.is_logged_in():
            raise CredError("这份登录状态已经过期了，请在桌面版里重新登录一次、"
                            "再把它导出的 login-state.json 拿来。")
        say("登录状态有效。")
        return session, profile
    except BaseException:
        try:
            session.close()
        except Exception:  # noqa: BLE001
            pass
        shutil.rmtree(profile, ignore_errors=True)
        raise


def login(username: str, password: str, log=None) -> tuple:
    """用用户自己的 iGEM 账号密码登录，返回一个已登录的 Session。"""
    say = log or (lambda _m: None)
    profile = Path(tempfile.mkdtemp(prefix="igem-web-"))
    session = igem.Session(headless=HEADLESS, profile_dir=profile, verbose=False,
                           timeout=240.0)
    session.start()
    try:
        page = session.page
        say("打开 iGEM 登录页…")
        page.goto(SIGN_IN, wait_until="domcontentloaded", timeout=60000)
        page.wait_for_timeout(900)

        if "accounts.igem.org" in page.url:
            try:
                page.fill("#identifier", username, timeout=25000)
                page.fill("#password", password, timeout=25000)
            except Exception as exc:  # noqa: BLE001
                raise CredError("找不到登录框，iGEM 登录页可能改版了。") from exc

            say(f"提交账号：{username}（{len(username)} 个字符），"
                f"密码 {len(password)} 个字符")
            say("填写完成，提交登录…")
            try:
                page.click("button[type=submit]", timeout=25000)
            except Exception:  # noqa: BLE001
                page.press("#password", "Enter")      # 兜底：回车提交

            # 等它离开登录页；期间盯着页面上的报错
            deadline = 90.0
            waited = 0.0
            alert = ""
            while waited < deadline:
                if "accounts.igem.org" not in page.url:
                    break
                page.wait_for_timeout(500)
                waited += 0.5
                if not alert:
                    alert = _visible_alert(page)
                    if alert:
                        hint = ""
                        low = alert.lower()
                        for key, text in _ALERT_HINTS.items():
                            if key in low:
                                hint = "\n" + text
                                break
                        raise CredError(f"iGEM 登录页提示：{alert}{hint}")
                if waited == 5.0:
                    say("还在登录页，等一下…")

            if "accounts.igem.org" in page.url:
                hint = _page_hint(page)
                raise CredError(
                    "登录没有成功。" + (f"（{hint}）" if hint else "")
                    + f"\n停留页面：{page.url}")

            say(f"已跳出登录页：{page.url}")

        if "accounts.igem.org" in page.url:
            raise CredError("登录后仍停在 iGEM 账号页，没有回到 registry。")

        say("核对登录状态…")
        if not session.is_logged_in():
            hint = _page_hint(page)
            raise CredError("登录提交成功但没检测到登录态。"
                            + (f"（{hint}）" if hint else "")
                            + f"\n当前页面：{page.url}")

        return session, profile          # profile 交给调用方删
    except BaseException:
        try:
            session.close()
        except Exception:  # noqa: BLE001
            pass
        shutil.rmtree(profile, ignore_errors=True)
        raise


# --------------------------------------------------------------------------
# 解析 / 分组 / 转换
# --------------------------------------------------------------------------

def _decode_dna(data_b64: str):
    """上传的 base64 → 临时文件 → 解析出 (parsed, pieces)。"""
    try:
        raw = base64.b64decode(data_b64 or "", validate=False)
    except Exception as exc:  # noqa: BLE001
        raise ValueError("上传的文件读不出来。") from exc
    if not raw:
        raise ValueError("上传的文件是空的。")

    fd, tmp_name = tempfile.mkstemp(suffix=".dna")
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "wb") as fh:      # 关掉句柄，否则 Windows 删不掉
            fh.write(raw)
        return engine.plan(tmp, require_features=True, rev_comp=True, log=None)
    finally:
        tmp.unlink(missing_ok=True)


def _piece_label(p) -> str:
    """分段列表里显示的名字（不带正向/反向后缀）。"""
    return p.title.rsplit("（", 1)[0]


def plan_segments(body: dict) -> list:
    """只解析、不登录：返回自动分好的段，给页面做「可合并」的列表。"""
    parsed, pieces = _decode_dna(body.get("data") or "")
    out = []
    for i, p in enumerate(pieces):
        out.append({
            "i": i,
            "label": _piece_label(p),
            "strand": p.strand,
            "coords": p.coords,
            "start": p.ranges[0][0] + 1,      # 显示用 1-based
            "end": p.ranges[-1][1],
            "length": len(p.seq),
            "kind": p.kind,
        })
    return out


def _rows_from_groups(parsed, pieces, groups):
    """按页面的分组把片段拼成输出行；groups 为 None 时每段一行。"""
    if groups is None:
        return [{"title": p.title, "seq": p.seq} for p in pieces]

    seen = set()
    rows = []
    for g in groups:
        if not isinstance(g, list) or not g:
            raise ValueError("分组格式不对。")
        idxs = []
        for v in g:
            try:
                i = int(v)
            except (TypeError, ValueError):
                raise ValueError("分组格式不对。")
            if i < 0 or i >= len(pieces) or i in seen:
                raise ValueError("分组里有重复或越界的片段。")
            seen.add(i)
            idxs.append(i)

        group = [pieces[i] for i in idxs]
        if len(group) == 1:
            p = group[0]
            rows.append({"title": p.title, "seq": p.seq})
            continue

        # 合并：按基因组坐标顺序把正链序列连起来
        start = group[0].ranges[0][0]
        end = group[-1].ranges[-1][1]
        name = " + ".join(_piece_label(p) for p in group)
        if len(name) > 48:
            name = f"合并({start + 1}..{end})"
        seq = "".join(parsed.seq[a:b] for p in group for (a, b) in p.ranges)
        rows.append({"title": f"{name}（正向）", "seq": seq})

    if len(seen) != len(pieces):
        raise ValueError("分组没有覆盖全部片段，会漏序列。")
    return rows


def convert(body: dict, emit) -> None:
    """把一次上传转成 txt，过程中用 emit() 往回调进度。"""
    username = (body.get("username") or "").strip()
    password = body.get("password") or ""
    filename = body.get("filename") or "input.dna"
    state_b64 = body.get("state") or ""

    # 两种登录方式二选一：账号密码，或者桌面版导出的 login-state.json
    if not state_b64 and (not username or not password):
        raise CredError("请填写 iGEM 账号和密码（或者上传 login-state.json）。")
    if state_b64:
        try:
            state_bytes = base64.b64decode(state_b64, validate=False)
        except Exception as exc:  # noqa: BLE001
            raise CredError("登录状态文件读不出来。") from exc

    def log(text: str) -> None:
        emit({"type": "log", "text": text})

    def step(index: int, total: int) -> None:
        emit({"type": "log", "text": f"查砖号 {index}/{total} …",
              "step": index, "total": total})

    # 先解析：不用登录就能知道文件有没有标注、分了几段
    parsed, pieces = _decode_dna(body.get("data") or "")
    rows = _rows_from_groups(parsed, pieces, body.get("groups"))

    emit({"type": "log", "text": f"{filename}：{len(parsed.seq)} bp，"
                                f"{len(rows)} 段，开始登录 iGEM…"})

    session = None
    profile = None
    try:
        if state_b64:
            session, profile = login_with_state(state_bytes, log=log)
        else:
            session, profile = login(username, password, log=log)
        emit({"type": "log", "text": "登录成功，开始查砖号…"})

        finder = engine.BrickFinder(profile_dir=profile, headless=True, verbose=False)
        finder.session = session            # 直接把已登录的会话塞进去

        found = finder.query_many([r["seq"] for r in rows], progress=log, on_step=step)

        lines = []
        total = len(rows)
        for i, r in enumerate(rows, 1):
            brick = found.get(r["seq"])
            body_text = f"{brick} {r['seq']}" if brick else r["seq"]
            lines.append(f"{r['title']}：{body_text}")
            emit({"type": "log",
                  "text": f"[{i}/{total}] {r['title']} → {brick or '查不到砖号，写序列'}"})

        emit({"type": "done", "lines": lines})
    finally:
        if session is not None:
            try:
                session.close()
            except Exception:  # noqa: BLE001
                pass
        if profile is not None:             # 真正用过的那个临时配置目录
            shutil.rmtree(profile, ignore_errors=True)


# --------------------------------------------------------------------------
# HTTP
# --------------------------------------------------------------------------

class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.0"

    def log_message(self, fmt, *args):      # 安静一点
        if VERBOSE:
            sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    # -- 小工具 ---------------------------------------------------------

    def _cors(self):
        self.send_header("Access-Control-Allow-Origin", ALLOW_ORIGIN)
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Methods", "POST, GET, OPTIONS")

    def _send(self, code, body: bytes, ctype: str):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self._cors()
        self.end_headers()
        self.wfile.write(body)

    def _json(self, code, obj):
        self._send(code, json.dumps(obj, ensure_ascii=False).encode("utf-8"),
                   "application/json; charset=utf-8")

    # -- 路由 -----------------------------------------------------------

    def do_OPTIONS(self):
        self.send_response(204)
        self._cors()
        self.end_headers()

    def do_GET(self):
        path = self.path.split("?")[0]
        if path in ("/", "/index.html"):
            page = HERE / "index.html"
            if not page.exists():
                self._json(500, {"error": "index.html 不见了"})
                return
            self._send(200, page.read_bytes(), "text/html; charset=utf-8")
        elif path == "/api/health":
            self._json(200, {"ok": True, "service": "igem-parts"})
        else:
            self._json(404, {"error": "not found"})

    def do_POST(self):
        path = self.path.split("?")[0]
        if path not in ("/api/convert", "/api/plan"):
            self._json(404, {"error": "not found"})
            return

        try:
            length = int(self.headers.get("Content-Length") or 0)
            raw = self.rfile.read(length) if length else b""
            body = json.loads(raw.decode("utf-8") or "{}")
        except Exception:  # noqa: BLE001
            self._json(400, {"error": "请求格式不对"})
            return

        # /api/plan：只解析，不登录，返回自动分好的段
        if path == "/api/plan":
            try:
                segments = plan_segments(body)
                self._json(200, {"segments": segments})
            except Exception as exc:  # noqa: BLE001
                if VERBOSE:
                    traceback.print_exc()
                self._json(200, {"error": f"{type(exc).__name__}: {exc}"})
            return

        # 流式返回：一行一个 JSON 事件
        self.send_response(200)
        self.send_header("Content-Type", "application/x-ndjson; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self._cors()
        self.end_headers()

        def emit(obj):
            try:
                self.wfile.write((json.dumps(obj, ensure_ascii=False) + "\n").encode("utf-8"))
                self.wfile.flush()
            except Exception:          # 客户端断了就算了
                pass

        with _RUN_LOCK:
            try:
                convert(body, emit)
            except CredError as exc:
                emit({"type": "error", "message": str(exc)})
            except Exception as exc:  # noqa: BLE001
                if VERBOSE:
                    traceback.print_exc()
                emit({"type": "error",
                      "message": f"{type(exc).__name__}: {exc}"})


def main() -> int:
    if sys.stdout is None:
        sys.stdout = open(os.devnull, "w", encoding="utf-8")
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass

    server = ThreadingHTTPServer((HOST, PORT), Handler)
    print(f"iGEM 砖号查询 · 网页版已启动： http://127.0.0.1:{PORT}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n已停止。", flush=True)
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
