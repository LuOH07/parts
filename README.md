# iGEM-parts

把 SnapGene 的 `.dna` 质粒文件转成「一行一段」的 iGEM 砖号清单。

- 左边上传 `.dna`、填自己的 iGEM 账号密码、点「开始查询」
- 右边看结果，一个按钮下载 `.txt`
- 全部功能在同一个页面里，不跳转

输出格式（一行一段，按质粒坐标从 1 号碱基排到最后一个）：

```
J23100（正向）：BBa_K823004
无标注区(96..112)（正向）：CGTACGTAGCTAGCTAG
B0034（反向）：AAAGAGGAGAAA
```

- 查得到砖号就写砖号，查不到就写这一段的序列 —— 所以整条质粒的每一段都在里面，不会漏
- 只认 100% 相同、而且长度完全一致的砖
- `type=source` 的标注、以及把别的标注整段包住的「区域」标注会跳过（它们不是部件）
- 引物结合位点不算砖

---

## ⚠️ 为什么必须有这个后端

iGEM 的接口只允许 `https://registry.igem.org` 这一个来源访问。实测：

```
GET https://api.registry.igem.org/v1/accounts/me
Origin: https://registry.igem.org
→ 200  access-control-allow-origin: https://registry.igem.org
        access-control-allow-credentials: true

Origin: https://luoh07.github.io
→ （不返回 access-control-allow-origin）
```

浏览器会因此拦掉请求。所以**单纯放一个静态页面在 GitHub Pages 上是没法登录和查 BLAST 的**，
登录和查询这两步只能放在服务端做。这个仓库里的 `app.py` 就是那个服务端。

好消息是：同一个 `app.py` 既当后端，也顺便把 `index.html` 发出去，
所以本地跑起来之后，页面和后端同源，什么都不用配。

---

## 本地运行（最简单）

需要 Python 3.10+，第一次装依赖：

```bash
pip install -r requirements.txt
python -m playwright install chromium
```

然后：

```bash
python app.py
```

打开 <http://127.0.0.1:8765> 就能用了。

Windows 上直接双击 **`启动网页版.bat`**（它会先启动服务，再自动打开浏览器）。

---

## 部署成公网链接

`app.py` 是标准 WSGI/HTTP 服务，读环境变量 `PORT`，可以直接丢到任何支持 Python 的平台。

**Render（免费档）**

1. 新建 Web Service，连这个仓库
2. Build Command：`pip install -r requirements.txt && python -m playwright install --with-deps chromium`
3. Start Command：`python app.py`
4. 部署完拿到的地址就是页面地址（页面和后端同源）

**Railway / Fly / 自己的服务器** 同理。

**只想把页面放 GitHub Pages？**
那就把 `index.html` 里这一行改成后端地址，页面放在 Pages、后端在别处：

```js
const API_BASE = "https://你的后端地址";
```

不填（留空）就表示页面和后端同源。

---

## 关于账号密码

- 每个使用者填**自己的** iGEM 账号密码，不是别人的
- 密码只在这一次请求里用于登录，不写日志、不落盘、服务端不保存
- 勾选「记住账号和密码」后，只保存在**使用者自己浏览器的 localStorage** 里，不会上传
- 每个请求各开一个临时浏览器配置目录，用完就删，请求之间不共享登录态

---

## 目录

```
index.html              单页界面（左边输入、右边输出、下载按钮）
app.py                  后端：发页面 + /api/convert（登录 iGEM、查 BLAST、流式返回）
requirements.txt
启动网页版.bat           Windows 一键启动
app/
  dna_to_bricks.py      .dna 解析、分段、输出格式
  igem_brick_finder.py  用 Playwright 操作 registry.igem.org 的 BLAST
```

## 限制

- 一个质粒一段一段查，十几段大概几分钟；页面右上角会实时显示进度
- 很短的段（6~12 bp 的 RBS / scar）BLAST 服务器不返回结果，只能写序列
- 一个进程同一时间只跑一个转换（Playwright 不适合并发抢同一个浏览器）
