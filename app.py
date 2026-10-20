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

# 同一个进程里同时跑多个 Playwright 会打架，串起来执行
_RUN_LOCK = threading.Lock()


class CredError(RuntimeError):
    """账号密码不对 / 登录没成功。"""


# --------------------------------------------------------------------------
# 登录
# --------------------------------------------------------------------------

def login(username: str, password: str) -> "igem.Session":
    """用用户自己的 iGEM 账号密码登录，返回一个已登录的 Session。"""
    profile = Path(tempfile.mkdtemp(prefix="igem-web-"))
    session = igem.Session(headless=True, profile_dir=profile, verbose=False,
                           timeout=240.0)
    session.start()
    try:
        page = session.page
        page.goto(SIGN_IN, wait_until="domcontentloaded", timeout=60000)
        page.wait_for_timeout(900)

        if "accounts.igem.org" in page.url:
            try:
                page.fill("#identifier", username, timeout=25000)
                page.fill("#password", password, timeout=25000)
                page.click("button[type=submit]", timeout=25000)
            except Exception as exc:  # noqa: BLE001
                raise CredError("找不到登录框，iGEM 登录页可能改版了。") from exc

            deadline = 240.0
            waited = 0.0
            while waited < deadline:
                if "accounts.igem.org" not in page.url:
                    break
                page.wait_for_timeout(500)
                waited += 0.5
            else:
                raise CredError("登录没有成功：账号或密码不对，或者这个账号需要额外验证"
                                "（两步验证 / 未激活）。")

        if "accounts.igem.org" in page.url or not session.is_logged_in():
            raise CredError("登录没有成功：账号或密码不对，或者这个账号需要额外验证。")

        return session
    except BaseException:
        try:
            session.close()
        except Exception:  # noqa: BLE001
            pass
        shutil.rmtree(profile, ignore_errors=True)
        raise


# --------------------------------------------------------------------------
# 转换
# --------------------------------------------------------------------------

def convert(body: dict, emit) -> None:
    """把一次上传转成 txt，过程中用 emit() 往回调进度。"""
    username = (body.get("username") or "").strip()
    password = body.get("password") or ""
    filename = body.get("filename") or "input.dna"
    data_b64 = body.get("data") or ""

    if not username or not password:
        raise CredError("请填写 iGEM 账号和密码。")
    try:
        raw = base64.b64decode(data_b64, validate=False)
    except Exception as exc:  # noqa: BLE001
        raise ValueError("上传的文件读不出来。") from exc
    if not raw:
        raise ValueError("上传的文件是空的。")

    def log(text: str) -> None:
        emit({"type": "log", "text": text})

    def step(index: int, total: int) -> None:
        emit({"type": "log", "text": f"查砖号 {index}/{total} …",
              "step": index, "total": total})

    # 先本地解析：不用登录就能知道文件有没有标注
    fd, tmp_name = tempfile.mkstemp(suffix=".dna")
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "wb") as fh:      # 关掉句柄，否则 Windows 删不掉
            fh.write(raw)
        parsed, pieces = engine.plan(tmp, require_features=True, rev_comp=True,
                                     log=log)
    finally:
        tmp.unlink(missing_ok=True)

    emit({"type": "log", "text": f"{filename}：{len(parsed.seq)} bp，"
                                f"{len(pieces)} 段，开始登录 iGEM…"})

    profile = Path(tempfile.mkdtemp(prefix="igem-web-"))
    session = None
    try:
        session = login(username, password)
        session.profile_dir = profile
        emit({"type": "log", "text": "登录成功，开始查砖号…"})

        finder = engine.BrickFinder(profile_dir=profile, headless=True, verbose=False)
        finder.session = session            # 直接把已登录的会话塞进去

        found = finder.query_many([p.seq for p in pieces], progress=log, on_step=step)
        engine.apply_bricks(pieces, found, log=log)
        lines = engine.output_lines(pieces)
        emit({"type": "done", "lines": lines})
    finally:
        if session is not None:
            try:
                session.close()
            except Exception:  # noqa: BLE001
                pass
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
        if self.path.split("?")[0] != "/api/convert":
            self._json(404, {"error": "not found"})
            return

        try:
            length = int(self.headers.get("Content-Length") or 0)
            raw = self.rfile.read(length) if length else b""
            body = json.loads(raw.decode("utf-8") or "{}")
        except Exception:  # noqa: BLE001
            self._json(400, {"error": "请求格式不对"})
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
