# TFBS Visualization Web Service

基于 `visualize_tfbs.R` 的 Web 封装，提供浏览器端的 TFBS / PlantCARE 可视化服务。

## 快速开始

### Docker 部署（推荐）

```bash
# 从 visualize_tfbs/ 根目录构建
docker build -t tfbs-web -f web/Dockerfile .

# 启动
docker run -d -p 6820:6820 --name tfbs-web tfbs-web

# 或使用 docker-compose
cd web
docker-compose up -d
```

访问 `http://localhost:6820`。

### 本地开发

```bash
cd web
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

# macOS / Linux
R_SCRIPT=../visualize_tfbs.R .venv/bin/uvicorn app:app --reload --port 6820
```

## 架构

```
浏览器 (HTML + vanilla JS + Tailwind CDN)
  │  POST /api/visualize (multipart form)
  ▼
FastAPI (web/app.py)
  │  subprocess → Rscript visualize_tfbs.R
  ▼
R 脚本输出 (PNG / PDF / JSON / CSV)
  │
  ▼
JSON 响应 { job_id, png_url, summary, downloads }
  │
  ▼
浏览器 fetch png_url → 显示图片 + 下载链接
```

**无缓存**：每个请求创建独立临时目录，180 秒后自动清理。

## API

### `POST /api/visualize`

提交文件并运行分析。

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| mode | string | 是 | `"tfbs"` 或 `"plantcare"` |
| fasta_file | file | TFBS 必需 | FASTA 序列文件 |
| prediction_file | file | TFBS 必需 | 预测结果 TSV |
| plantcare_file | file | PlantCARE 必需 | PlantCARE tab 文件 |
| single_strand | bool | 否 | 合并正负链，默认 false |
| pvalue_threshold | float | TFBS 可选 | p 值阈值，默认 1e-6 |
| motif_filter | string | PlantCARE 可选 | 逗号分隔的 motif 白名单 |

响应示例：

```json
{
  "job_id": "a1b2c3d4e5f6",
  "png_url": "/api/result/a1b2c3d4e5f6/output.png",
  "summary": { ... },
  "downloads": {
    "png": "/api/result/a1b2c3d4e5f6/output.png",
    "pdf": "/api/result/a1b2c3d4e5f6/output.pdf",
    "coordinates": "/api/result/a1b2c3d4e5f6/coordinates.csv",
    "summary_json": "/api/result/a1b2c3d4e5f6/summary.json",
    "bundle": "/api/result/a1b2c3d4e5f6/bundle.zip"
  }
}
```

PlantCARE 模式额外提供 `raw_mapping` 下载链接。

### `GET /api/result/{job_id}/{filename}`

下载单个输出文件。有效的文件名：`output.png`、`output.pdf`、`coordinates.csv`、`summary.json`、`raw_mapping.csv`。

### `GET /api/result/{job_id}/bundle.zip`

下载全部结果的 ZIP 包。

## 项目结构

```
web/
  Dockerfile              # 单阶段构建（Python + R）
  docker-compose.yml      # 端口映射 6820
  .dockerignore
  requirements.txt        # fastapi, uvicorn, python-multipart
  app.py                  # FastAPI 后端
  static/
    index.html            # 前端页面（Tailwind CDN）
    app.js                # 前端逻辑（vanilla JS）
    style.css             # 补充样式
```

## 设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 前端框架 | 纯 HTML + vanilla JS | 无需构建步骤，Tailwind 通过 CDN 引入 |
| 响应格式 | JSON（PNG 通过 URL 二次获取） | 避免 header 编码问题，结构清晰 |
| Job 存储 | 内存 dict + /tmp 临时目录 | 无外部依赖，10 分钟 TTL 足够 |
| 清理机制 | asyncio 后台任务 | 进程内完成，无需 cron |
| R 调用方式 | 同步 subprocess | R 脚本通常 <5s，无需任务队列 |
| Docker 构建 | 单阶段 | 无前端构建步骤，不需要 Node.js |
| 配置 | 始终 `--no-config` | 容器内无 config 文件 |
| 文件大小限制 | 10MB/文件 | 启动子文件极小，10MB 已很宽松 |
| 超时 | 120 秒 | 防止僵尸进程 |
