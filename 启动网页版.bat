@echo off
chcp 65001 >nul
cd /d "%~dp0"

REM  iGEM brick lookup - web version
REM  starts the local server and opens the page in your browser

set PY=
if exist "python\python.exe" set PY=python\python.exe
if not defined PY if exist "..\iGEM砖号查询\python\python.exe" set PY=..\iGEM砖号查询\python\python.exe
if not defined PY set PY=python

echo 启动中… 浏览器会自动打开 http://127.0.0.1:8765
start "iGEM brick lookup server" "%PY%" app.py
timeout /t 3 >nul
start "" http://127.0.0.1:8765
