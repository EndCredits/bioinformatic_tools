# visualize_tfbs.R

启动子区域转录因子结合位点 (TFBS) / 顺式作用元件可视化工具。支持两种输入格式：PlantRegMap TFBS 预测结果和 PlantCARE 顺式作用元件预测结果，输出 publication-ready 的 PDF/PNG 图谱和 JSON 统计摘要。

## 引用

本项目使用 PlantRegMap 和 PlantCARE 的成果，如果您使用了上述任意一种方式预测，请记得规范引用:

```
PlantRegMap:
Tian, F., Yang, D.C., Meng, Y.Q., Jin, J. and Gao, G. (2020) PlantRegMap: charting functional regulatory maps in plants. Nucleic Acids Res 48, D1104-D1113.
Jin, J., Tian, F., Yang, D.C., Meng, Y.Q., Kong, L., Luo, J. and Gao, G. (2017) PlantTFDB 4.0: toward a central hub for transcription factors and regulatory interactions in plants. Nucleic Acids Res, 45, D1040-D1045.
Jin, J., He, K., Tang, X., Li, Z., Lv, L., Zhao, Y., Luo, J. and Gao, G. (2015) An Arabidopsis transcriptional regulatory map reveals distinct functional and evolutionary features of novel transcription factors. Mol Biol Evol, 32, 1767-1773.

PlantCARE

Lescot, M., Déhais, P., Thijs, G., Marchal, K., Moreau, Y., Van de Peer, Y., Rouzé, P., & Rombauts, S. (2002). PlantCARE, a database of plant cis-acting regulatory elements and a portal to tools for in silico analysis of promoter sequences. Nucleic Acids Research, 30(1), 325-327. doi:10.1093/nar/30.1.325
```

## 依赖

**R 包**（必需）：

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "viridis", "scales", "jsonlite"))
```

**R 包**（可选，仅 `--config` 需要）：

```r
install.packages("yaml")
```

## 快速开始

```bash
# TFBS 模式（PlantRegMap / FIMO 输出）
Rscript visualize_tfbs.R promoter.fasta predictions.txt output_prefix

# PlantCARE 模式
Rscript visualize_tfbs.R --plantcare plantCARE_output.tab promoter.fa output_prefix

# 单链视图（正负链合并显示）
Rscript visualize_tfbs.R --single promoter.fasta predictions.txt output_prefix
```

## 输入格式

### TFBS 模式（默认）

需要两个输入文件：

**FASTA 文件** — 启动子序列：

```
>gene_name
ATCGATCGATCG...
```

**预测文件** — Tab 分隔，首行为 `#` 开头的列名行，必需列：

```
#pattern name  Family  sequence name  start  stop  strand  score  p-value  q-value  matched sequence
Aradu.D4D3B    TCP     UGT3A-promoter  542    571   -       23.92  7.3e-09  2.6e-05  GGTGGGCCC...
```

列名不区分大小写，`-` 和 `.` 会自动转为 `_`。必需列：`family`, `start`, `stop`, `strand`, `score`, `p_value`, `q_value`。

### PlantCARE 模式 (`--plantcare`)

**PlantCARE 文件** — Tab 分隔，每行一个元件，无表头（表头行会自动检测并跳过）：

```
seq_name  family  motif_sequence  start  length  strand  organism  description
Gene1     ABRE3a  TACGTG          316    6       +       Zea mays  ABA responsiveness
```

共 7 列，第 8 列起为可选的 description 信息。

**FASTA 文件**（可选）— 传 `-` 跳过，脚本会从数据中推断序列长度。

## 命令行选项

| 选项 | 说明 |
|------|------|
| `--plantcare` | 使用 PlantCARE 输入格式 |
| `--single` | 正负链合并为单行视图（两种模式均支持） |
| `--exclude MOTIFS` | 逗号分隔的 motif 名称列表，覆盖配置文件中的排除列表 |
| `--no-exclude` | 禁用所有排除规则，显示全部 motif |
| `--config FILE` | 指定 YAML 配置文件路径 |
| `--no-config` | 忽略配置文件，仅使用内置默认值 |

### 过滤优先级

1. **排除列表** (`--exclude` 或配置文件) — 先执行
2. **Motif 白名单**（第 4 个位置参数） — 后执行
3. **单例过滤** — 移除仅出现一次的 motif 家族

`--exclude` 覆盖配置文件中的排除列表。白名单无匹配时回退到排除+单例过滤后的数据。

## 示例

```bash
# TFBS 模式，指定 p-value 阈值
Rscript visualize_tfbs.R promoter.fasta predictions.txt output 1e-6

# PlantCARE 模式，不需要 FASTA
Rscript visualize_tfbs.R --plantcare plantCARE_output.tab - output

# PlantCARE 模式，只显示特定 motif
Rscript visualize_tfbs.R --plantcare plantCARE_output.tab promoter.fa output TATA-box,ABRE,G-box

# 排除特定 motif
Rscript visualize_tfbs.R --plantcare --exclude "G-box,ABRE" plantCARE_output.tab promoter.fa output

# 单链视图（TFBS 或 PlantCARE 均可）
Rscript visualize_tfbs.R --single promoter.fasta predictions.txt output
Rscript visualize_tfbs.R --plantcare --single plantCARE_output.tab promoter.fa output
```

## 输出文件

| 文件 | 内容 |
|------|------|
| `<prefix>.pdf` | 主图（矢量，适合论文） |
| `<prefix>.png` | PNG 版本（300 dpi） |
| `<prefix>_summary.json` | 统计摘要（序列长度、元件数量、分布、富集区域等） |
| `<prefix>_coordinates.csv` | 坐标变换后的元件位置 |
| `<prefix>_raw_mapping.csv` | 原始 PlantCARE 数据（仅 PlantCARE 模式） |

### strand 信息说明

- **双链模式**：`_coordinates.csv` 和 `_summary.json` 中 strand 字段为原始 `+` / `-`
- **单链模式** (`--single`)：`_coordinates.csv` 保留原始 strand，`_summary.json` 中 `enriched_regions` 的 strand 为聚合后的值（如 `+/-` 表示正负链元件在同一区域重叠）

## 显著性标记

### TFBS 模式

基于 p-value：

| 标记 | 阈值 |
|------|------|
| `***` | p < 1e-8 |
| `**` | p < 1e-7 |
| `*` | p < 1e-6 |

使用 CLI 指定 p-value 阈值时，标记阈值自动对齐（`low` = 用户阈值，`medium` = low/10，`high` = low/100）。

### PlantCARE 模式

基于聚合后元件的出现次数（count）：

| 标记 | 阈值 |
|------|------|
| `***` | count >= enriched_min_count × 3 |
| `**` | count >= enriched_min_count × 2 |
| `*` | count >= enriched_min_count |

## 配置文件

脚本自动加载同目录下的 `tfbs_config.yaml`。配置文件格式：

```yaml
# 排除的 motif 列表
exclude_motifs:
  - "Unnamed__1"
  - "TATA-box"
  - "CAAT-box"

# 移除仅出现一次的 motif 家族
remove_singletons: true

# 显著性阈值（TFBS 模式）
significance:
  high: 1e-8     # ***
  medium: 1e-7   # **
  low: 1e-6      # *

# "n=x" 标签的最小出现次数
enriched_min_count: 3

# 可视化参数
strand_y_offset: 0.4      # 双链模式 DNA 链的 Y 轴间距
stripe_width: 200         # 背景条纹宽度 (bp)
label_min_dist_ratio: 0.045  # 标签最小间距（序列长度的比例）
alpha_score_range: [4, 12]   # PlantCARE alpha 映射范围（motif 长度）
lane_height: 0.15         # 标签泳道高度
merge_threshold: 0        # 位点聚合间距阈值 (bp)
```

## 透明度

图中每个位点的透明度反映其显著性/长度：

- **TFBS 模式**：p-value 越小 → 越不透明（alpha 0.3–1.0）
- **PlantCARE 模式**：motif 越长 → 越不透明（自动适应数据范围）

## 图谱元素

- **红色虚线** — TSS 位置（位置 1）
- **黑色实线** — DNA 链（双链模式两条，单链模式一条）
- **彩色点** — 预测位点，颜色按 TF 家族/motif 区分
- **彩色标签** — 家族名称，自动防重叠分配泳道
- **红色标记** — `***` / `**` / `*` 显著性
- **灰色斜体** — `n=x` 富集区域标记（count >= enriched_min_count）
- **背景条纹** — 200 bp 间隔，便于定位
