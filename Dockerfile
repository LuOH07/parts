FROM python:3.11-slim

ENV PYTHONUNBUFFERED=1 \
    PLAYWRIGHT_BROWSERS_PATH=/ms-playwright \
    HEADLESS=1

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt \
    && python -m playwright install --with-deps chromium \
    && apt-get update \
    && apt-get install -y --no-install-recommends xvfb \
    && rm -rf /var/lib/apt/lists/*

COPY . .

ENV PORT=8765
EXPOSE 8765

# 默认 headless（512MB 免费档撑不住有头 Chromium）。
# 想试有头：在平台的环境变量里把 HEADLESS 改成 0 即可 —— xvfb 已经装好，
# 下面这条 xvfb-run 就会给它一块虚拟屏幕，不用重新构建。
CMD ["xvfb-run", "-a", "-s", "-screen 0 1366x900x24", "python", "app.py"]
