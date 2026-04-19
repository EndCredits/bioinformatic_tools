# TFBS Visualization Web Service

基于 `visualize_tfbs.R` 的 Web 封装，提供浏览器端的 TFBS / PlantCARE 可视化服务。

## 快速开始

### Nix（推荐）

```bash
# 开发环境
cd visualize_tfbs/web
nix develop

# 启动服务
R_SCRIPT=../visualize_tfbs.R uvicorn app:app --reload --port 6820
```

### Docker

```bash
# 从 visualize_tfbs/ 根目录构建
docker build -t tfbs-web -f web/Dockerfile .
docker run -d -p 6820:6820 tfbs-web

# 或使用 docker-compose
cd web && docker-compose up -d
```

### Nix OCI 镜像（无需 Docker）

```bash
cd visualize_tfbs/web
nix build .#default          # 构建 OCI tarball
docker load < result         # 加载到 Docker
# 或
skopeo copy docker-archive:result docker://registry/tfbs-web:latest
```

访问 `http://localhost:6820`。

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

**无缓存**：每个请求创建独立临时目录，R 脚本重新生成输出，180 秒后自动清理。

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
  flake.nix              # Nix dev shell + OCI 镜像构建
  flake.lock             # 依赖锁定（nixpkgs commit）
  Dockerfile             # Docker 构建（备选方案）
  docker-compose.yml     # 端口映射 6820
  .dockerignore
  requirements.txt       # Python 依赖声明
  requirements.lock      # Python 依赖锁定（pip freeze）
  r-packages.lock        # R 依赖锁定
  app.py                 # FastAPI 后端
  static/
    index.html           # 前端页面（Tailwind CDN）
    app.js               # 前端逻辑（vanilla JS）
    style.css            # 补充样式
```

## 依赖管理

项目支持两种依赖管理方式：

### Nix（精确可重现）

由 `flake.lock` 中的 nixpkgs commit 哈希精确锁定。所有子 flake 的 nixpkgs 版本通过脚本同步：

```bash
./scripts/nix-flake-update.sh --check   # 检查一致性
./scripts/nix-flake-update.sh           # 全部同步
```

### Docker / pip（版本锁定）

- `requirements.lock` — Python 包版本（`pip freeze`）
- `r-packages.lock` — R 包版本

## 设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 前端框架 | 纯 HTML + vanilla JS | 无需构建步骤，Tailwind 通过 CDN 引入 |
| 响应格式 | JSON（PNG 通过 URL 二次获取） | 避免 header 编码问题，结构清晰 |
| Job 存储 | 内存 dict + /tmp 临时目录 | 无外部依赖，180 秒 TTL 足够 |
| 清理机制 | asyncio 后台任务 | 进程内完成，无需 cron |
| R 调用方式 | 同步 subprocess | R 脚本通常 <5s，无需任务队列 |
| 依赖管理 | Nix flake + Docker 备选 | Nix 精确可重现，Docker 作为兼容方案 |
| 配置 | 始终 `--no-config` | 容器内无 config 文件 |
| PNG 后端 | Cairo（`type = "cairo"`） | 跨平台兼容，避免 Nix 下 bus error |
| 文件大小限制 | 10MB/文件 | 启动子文件极小，10MB 已很宽松 |
| 超时 | 120 秒 | 防止僵尸进程 |
