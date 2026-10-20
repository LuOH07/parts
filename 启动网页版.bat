@echo off
chcp 65001 >nul
cd /d "%~dp0"

REM  iGEM brick lookup - web version
REM  starts the local server, then opens the page in your browser

set "PY="
if exist "python\python.exe" set "PY=python\python.exe"
if not defined PY for /d %%d in ("..\iGEM*") do if exist "%%d\python\python.exe" set "PY=%%d\python\python.exe"
if not defined PY set "PY=python"

echo.
echo   iGEM brick lookup / web version
echo   python : %PY%
echo   page   : http://127.0.0.1:8765
echo.
echo   Starting the local server, your browser will open in a moment.
echo   Keep this window open while you use it. Close it to stop.
echo.

start "igem-web-server" "%PY%" app.py
timeout /t 3 >nul
start "" http://127.0.0.1:8765
