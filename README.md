# iGEM 砖号查询（网页版）

https://github.com/LuOH07/parts

把 SnapGene 的 `.dna` 质粒文件转成「一行一段」的 iGEM 砖号清单。

- 左边上传 `.dna`、填自己的 iGEM 账号密码、点「开始查询」
- 中间那栏会把 SnapGene 自动分好的段列出来，**勾相邻的几段可以手动合并**，
  合并后的整条序列会拿去查砖号（适合“一个部件被拆成了好几条标注”的情况）
- 右边看结果，一个按钮下载 `.txt`
- 全部功能在同一个页面里，不跳转

输出格式（一行一段，砖号后面跟这一段对应的序列）：

```
J23100（正向）：BBa_K823004 TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC
无标注区(96..112)（正向）：CGTACGTAGCTAGCTAG
B0034（反向）：AAAGAGGAGAAA
```

- 查得到砖号就写「砖号 + 序列」，查不到就只写这一段的序列 —— 每一段都在里面，不会漏
- 反向链的 feature 输出的是它自己的 5'→3'（反向互补），标签写「反向」
- 手动合并出来的段，按基因组坐标顺序取正链序列，标签写「正向」
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

## 一键部署（推荐，全程只用浏览器）

点下面的按钮，用 GitHub 登录一下，剩下的它自己做完：

- **Render**（推荐，带 `render.yaml` 蓝图）— <https://render.com/deploy?repo=https://github.com/LuOH07/parts>
- **Railway** — 打开 <https://railway.app/new> → Deploy from GitHub repo → 选 `LuOH07/parts`（会自动按 `Dockerfile` 构建）

仓库里已经带好 `Dockerfile`、`render.yaml`、`Procfile`，不用手填 Build / Start 命令。
部署完成给你的那个网址 **本身就是页面**（页面和后端同源），发给谁都能直接用。

免费档的注意事项：空转会休眠，第一次打开等十几秒；单个请求有平台超时，
段数很多的质粒建议分几次跑，或者自己升一档。

### 手动部署（任何支持 Python / Docker 的平台）

| 方式 | 命令 |
| --- | --- |
| Docker | `docker build -t igem-parts . && docker run -p 8765:8765 igem-parts` |
| 直接跑 | `pip install -r requirements.txt` + `python -m playwright install --with-deps chromium` + `python app.py` |

服务读环境变量 `PORT`，平台给什么就用什么。

### 只想把页面放 GitHub Pages？

把 `index.html` 里这一行改成后端地址，页面放 Pages、后端放别处：

```js
const API_BASE = "https://你的后端地址";
```

留空表示页面和后端同源。

---

## 在浏览器里用（部署完之后）

1. 打开后端给的那个网址
2. 左边拖入 `.dna`
3. 填自己的 iGEM 账号密码（可以勾「记住」）
4. 点「开始查询」，右边实时出进度，跑完点「下载 txt」

不需要在本机装任何东西。

---

## 本机直接跑（只是方便调试，不是必须）

```bash
pip install -r requirements.txt
python -m playwright install chromium
python app.py          # 打开 http://127.0.0.1:8765
```

Windows 上双击 `启动网页版.bat` 也行。

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
