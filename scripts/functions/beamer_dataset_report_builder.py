#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""为CRLM休眠课题自动构建模块化Beamer汇报源码。

脚本职责：
  1. 扫描当前数据集的结果目录，自动发现已有图片和CSV结果；
  2. 在temporary/beamer/generated_tables中生成Beamer专用表格预览；
  3. 在scripts/beamer/sections下生成分节tex文件；
  4. 生成scripts/beamer/beamer_report.tex主文件；
  5. 可选调用latexmk编译，并把最终PDF复制到results/reports/beamer。

设计原则：
  - 分析脚本负责产生真实结果，Beamer构建器只负责展示和排版；
  - 表格展示只做列筛选和预览，不改变原始CSV；
  - 图片引用优先使用PNG，避免在Beamer中插入大量矢量PDF导致编译过慢；
  - 展示顺序按“分析方案/交集方案”组织，便于汇报时沿同一证据链阅读。
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
from decimal import Decimal, InvalidOperation
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


# 项目根目录由当前脚本位置反推，避免在不同终端工作目录下运行时路径漂移。
PROJECT_ROOT = Path(__file__).resolve().parents[2]

# 当前构建的数据集。后续如果扩展到新数据集，只需抽象这两个常量和结果路径。
DATASET_ID = "GSE114012"
DATA_TYPE = "ngs"

RESULT_ROOT = PROJECT_ROOT / "results" / DATA_TYPE / DATASET_ID
PLOT_ROOT = RESULT_ROOT / "plots"
TABLE_ROOT = RESULT_ROOT / "tables"
INTERSECT_ROOT = RESULT_ROOT / "intersect"
TF_ROOT = RESULT_ROOT / "TF"
TF_SUMMARY_ROOT = RESULT_ROOT / "TF_summary"

BEAMER_ROOT = PROJECT_ROOT / "scripts" / "beamer"
SECTION_ROOT = BEAMER_ROOT / "sections"
PROJECT_SECTION_ROOT = SECTION_ROOT / "project"
DATASET_SECTION_ROOT = SECTION_ROOT / DATASET_ID
BUILD_ROOT = PROJECT_ROOT / "temporary" / "beamer"
GENERATED_TABLE_ROOT = BUILD_ROOT / "generated_tables" / DATASET_ID
REPORT_ROOT = PROJECT_ROOT / "results" / "reports" / "beamer"
MAIN_TEX = BEAMER_ROOT / "beamer_report.tex"
FINAL_PDF = REPORT_ROOT / "beamer_report.pdf"

THEME_DIR = PROJECT_ROOT / "data" / "templates" / "beamer" / "SimplePlus-BeamerTheme"


SCRIPT_TITLES = {
    "00": ("样本结构质控", "TPM 相关性确认 LRC/BULK 模型基础"),
    "01": ("LRC/BULK差异转录程序", "休眠样状态相对循环样状态的 DEG 证据"),
    "02": ("跨模型稳定候选基因", "从单模型 DEG 提炼可重复的休眠相关基因"),
    "04": ("多模型DEG方向总览", "并列比较不同细胞系的显著上/下调基因"),
    "05": ("Top DEG表达分离", "用样本层面热图验证候选基因集的状态区分能力"),
    "06": ("GSEA运行概览", "基于全量 ranked genes 的通路富集批量运算"),
    "07": ("GSEA图表配对解读", "同一富集设计的 dotplot 与统计表连续展示"),
    "08": ("TF富集方法框架", "六种 TF 富集/活性推断方法的输入与定位"),
    "09": ("TF候选整合与交集", "将多方法 TF 证据压缩为可验证候选"),
}

# Beamer中每个表格默认展示前21行；过长表格只预览核心列，完整CSV仍保留在结果目录。
TABLE_PREVIEW_ROWS = 21

# 生成大量表格预览时使用线程池并行，主要加速CSV读取和tex片段写入。
TABLE_PREVIEW_WORKERS = max(2, min((os.cpu_count() or 4), 12))
GENERATED_PREVIEW_TABLES: set[Path] = set()

ANALYSIS_ORDER_PRIORITY = {
    "ALL": 0,
    "DLD1_HCT15_SW48": 1,
    "DLD1_HCT15": 2,
    "DLD1": 10,
    "HCT15": 11,
    "HT55": 12,
    "RKO": 13,
    "SW48": 14,
    "SW948": 15,
}

RESULT_STEM_ORDER = {
    "summary": 0,
    "gene_list": 1,
    "significant_genes": 2,
    "deg_results": 3,
}

GSEA_OUTPUT_ORDER = {
    "hallmark": 0,
    "CP_BIOCARTA": 1,
    "CP_KEGG_MEDICUS": 2,
    "CP_KEGG_LEGACY": 3,
    "CP_REACTOME": 4,
    "CP_WIKIPATHWAYS": 5,
    "TFT_TFT_LEGACY": 6,
    "TFT_GTRD": 7,
    "GO_BP": 8,
    "GO_CC": 9,
    "GO_MF": 10,
    "HPO": 11,
    "C6": 12,
    "IMMUNESIGDB": 13,
}

GSEA_DISPLAY_LABELS = {
    "hallmark": "Hallmark",
    "CP_BIOCARTA": "BioCarta",
    "CP_KEGG_MEDICUS": "KEGG Medicus",
    "CP_KEGG_LEGACY": "KEGG Legacy",
    "CP_REACTOME": "Reactome",
    "CP_WIKIPATHWAYS": "WikiPathways",
    "TFT_TFT_LEGACY": "TFT Legacy",
    "TFT_GTRD": "TFT GTRD",
    "GO_BP": "GO Biological Process",
    "GO_CC": "GO Cellular Component",
    "GO_MF": "GO Molecular Function",
    "HPO": "Human Phenotype Ontology",
    "C6": "Oncogenic Signatures C6",
    "IMMUNESIGDB": "ImmuneSigDB",
}

TF_METHOD_ORDER = {
    "dorothea": 0,
    "chea3": 1,
    "viper": 2,
    "enrichr": 3,
    "trrust": 4,
    "collectri": 5,
}

TF_METHOD_LABELS = {
    "dorothea": "DoRothEA",
    "chea3": "ChEA3",
    "viper": "VIPER",
    "enrichr": "ENRICHR",
    "trrust": "TRRUST",
    "collectri": "CollecTRI",
}

TF_INTERSECTION_ORDER = {
    "ALL_6_METHODS": 0,
    "WITHOUT_CHEA3": 1,
    "DOROTHEA_CHEA3_VIPER": 2,
    "ENRICHR_TRRUST_COLLECTRI": 3,
    "ACTIVITY_METHODS_VIPER_COLLECTRI": 4,
    "ORA_METHODS_DOROTHEA_TRRUST": 5,
    "API_METHODS_CHEA3_ENRICHR": 6,
    "SIGNED_NETWORK_DOROTHEA_VIPER_COLLECTRI": 7,
    "LIST_BASED_ORA_API_DOROTHEA_CHEA3_ENRICHR_TRRUST": 8,
    "CHIP_LITERATURE_EVIDENCE_CHEA3_ENRICHR_TRRUST": 9,
    "BROAD_DATABASE_EVIDENCE_CHEA3_ENRICHR_COLLECTRI": 10,
    "CURATED_REGULON_NETWORK_DOROTHEA_TRRUST_COLLECTRI": 11,
}

TF_INTERSECTION_LABELS = {
    "ALL_6_METHODS": "六方法共同交集",
    "WITHOUT_CHEA3": "去除ChEA3后的五方法交集",
    "DOROTHEA_CHEA3_VIPER": "DoRothEA/ChEA3/VIPER交集",
    "ENRICHR_TRRUST_COLLECTRI": "ENRICHR/TRRUST/CollecTRI交集",
    "ACTIVITY_METHODS_VIPER_COLLECTRI": "活性推断方法交集",
    "ORA_METHODS_DOROTHEA_TRRUST": "ORA文献/调控网络交集",
    "API_METHODS_CHEA3_ENRICHR": "API数据库方法交集",
    "SIGNED_NETWORK_DOROTHEA_VIPER_COLLECTRI": "带方向调控网络交集",
    "LIST_BASED_ORA_API_DOROTHEA_CHEA3_ENRICHR_TRRUST": "列表富集/API方法综合交集",
    "CHIP_LITERATURE_EVIDENCE_CHEA3_ENRICHR_TRRUST": "ChIP/文献证据交集",
    "BROAD_DATABASE_EVIDENCE_CHEA3_ENRICHR_COLLECTRI": "广谱数据库证据交集",
    "CURATED_REGULON_NETWORK_DOROTHEA_TRRUST_COLLECTRI": "人工整理调控网络交集",
}

# Beamer表格预览列配置。
# 这里故意只挑选汇报时最有用、版式最稳的列；如后续某类结果需要增删列，
# 只改这一处即可，不需要改R脚本或分析输出。
TABLE_COLUMN_PRESETS = {
    "deg_summary": [
        "Dataset",
        "Data_Type",
        "Analysis_Name",
        "Contrast",
        "Samples_Used",
        "Up",
        "Down",
        "Total_Significant_Genes",
        "P_Value_Column",
        "P_Value_Cutoff",
        "LogFC_Cutoff",
    ],
    "deg_gene_table": [
        "Symbol",
        "Ensembl",
        "Entrez",
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val",
    ],
    "intersect_summary": [
        "Selected_Analyses",
        "Total_Intersected_Genes",
        "Common_Up",
        "Common_Down",
        "Mixed_Direction",
    ],
    "intersect_gene_list": [
        "Symbol",
        "Ensembl",
        "Entrez",
    ],
    "intersect_deg_results": [
        "Symbol",
        "Ensembl",
        "Entrez",
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val",
    ],
    "gsea_summary": [
        "Analysis_Name",
        "GeneSet_Name",
        "Ranked_Genes",
        "GSEA_Terms",
        "Positive_NES",
        "Negative_NES",
        "Single_Pathway_Plots",
    ],
    "gsea_result": [
        "ID",
        "NES",
        "pvalue",
        "p.adjust",
        "qvalue",
    ],
    "tf_run_summary": [
        "Input_Type",
        "TF_Analysis_Name",
        "Methods_Integrated",
        "Intersection_Schemes",
    ],
    "tf_intersection_summary": [
        "Input_Type",
        "TF_Analysis_Name",
        "Intersection_Name",
        "Required_Methods",
        "Number_Of_Methods",
        "Intersected_TF_Count",
        "Reported_Top_N",
    ],
    "tf_candidates": [
        "Consensus_Rank",
        "TF",
        "Source_Method_Count",
        "Source_Methods",
        "CheA3_Library_Count",
        "CheA3_Integrated_TopRank",
        "DoRothEA",
        "ChEA3",
        "VIPER",
        "ENRICHR",
        "TRRUST",
        "CollecTRI",
    ],
    "tf_method_final_summary": [
        "Input_Type",
        "TF_Analysis_Name",
        "Method",
        "Final_TF_Count",
    ],
    "tf_method_final": [
        "Rank",
        "TF",
        "Score",
        "P_Value",
        "Adjusted_P_Value",
        "NES",
        "Direction",
        "CheA3_Library_Count",
        "CheA3_Integrated_TopRank",
    ],
}


def rel(path: Path | str) -> str:
    path = Path(path)
    try:
        return path.resolve().relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def tex_escape(text: object) -> str:
    text = str(text)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def sanitize_name(text: str, default: str = "section") -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip())
    text = re.sub(r"_+", "_", text).strip("_")
    return text or default


def sort_analysis_name(name: str) -> tuple[int, str]:
    return (ANALYSIS_ORDER_PRIORITY.get(name, 100), name)


def sort_result_csv_by_stem(csv_path: Path) -> tuple[int, str]:
    flat_path = normalize_result_csv_path(csv_path)
    return (RESULT_STEM_ORDER.get(flat_path.stem, 50), flat_path.as_posix())


def sort_gsea_csv(csv_path: Path) -> tuple[int, str]:
    flat_path = normalize_result_csv_path(csv_path)
    try:
        parts = flat_path.relative_to(TABLE_ROOT).parts
        geneset = parts[2]
    except ValueError:
        geneset = flat_path.parent.name
    priority = GSEA_OUTPUT_ORDER.get(geneset, 100)
    return (priority, geneset)


def gsea_display_label(geneset: str) -> str:
    return GSEA_DISPLAY_LABELS.get(geneset, geneset)


def tf_method_label(method: str) -> str:
    return TF_METHOD_LABELS.get(method, method)


def tf_intersection_label(scheme: str) -> str:
    return TF_INTERSECTION_LABELS.get(scheme, scheme)


def sort_tf_csv(csv_path: Path) -> tuple[int, str, str]:
    flat_path = normalize_result_csv_path(csv_path)
    parts = flat_path.parts
    if "method_final" in parts and flat_path.name == "method_final_summary.csv":
        return (0, "", flat_path.as_posix())
    if "method_final" in parts:
        method = flat_path.parent.name
        return (1, f"{TF_METHOD_ORDER.get(method, 99):02d}_{method}", flat_path.as_posix())
    if "intersections" in parts and flat_path.name == "intersection_summary.csv":
        return (2, "", flat_path.as_posix())
    if "intersections" in parts and flat_path.name == "summary.csv":
        scheme = flat_path.parent.parent.name
        return (3, f"{TF_INTERSECTION_ORDER.get(scheme, 99):02d}_{scheme}", flat_path.as_posix())
    if "intersections" in parts and flat_path.name == "top10_tf_candidates.csv":
        scheme = flat_path.parent.parent.name
        return (4, f"{TF_INTERSECTION_ORDER.get(scheme, 99):02d}_{scheme}", flat_path.as_posix())
    return (10, "", flat_path.as_posix())


def sort_tf_group_key(input_type: str, analysis: str) -> tuple[int, tuple[int, str]]:
    type_priority = {"deg": 0, "DEG": 0, "intersect": 1, "INTERSECT": 1}
    return (type_priority.get(input_type, 9), sort_analysis_name(analysis))


def files(root: Path, pattern: str) -> list[Path]:
    if not root.exists():
        return []
    return sorted(root.rglob(pattern))


def dirs(root: Path) -> list[Path]:
    if not root.exists():
        return []
    return sorted([p for p in root.iterdir() if p.is_dir()])


def normalize_result_csv_path(csv_path: Path) -> Path:
    """Return the current flat CSV path for old or new result layouts."""
    csv_path = Path(csv_path)
    if csv_path.parent.name in {"csv", "md", "tex"}:
        return csv_path.parent.parent / f"{csv_path.stem}.csv"
    return csv_path


def legacy_result_csv_path(csv_path: Path) -> Path:
    """Return the historical output_dir/csv/result.csv path for compatibility."""
    flat_path = normalize_result_csv_path(csv_path)
    return flat_path.parent / "csv" / flat_path.name


def resolve_result_csv_path(csv_path: Path) -> Path:
    """Prefer flat CSV; fall back to historical csv/ directory if needed."""
    flat_path = normalize_result_csv_path(csv_path)
    legacy_path = legacy_result_csv_path(flat_path)
    if flat_path.exists():
        return flat_path
    if legacy_path.exists():
        return legacy_path
    return flat_path


def find_result_csv(output_dir: Path, stem: str) -> Path | None:
    """Find output_dir/stem.csv while remaining compatible with output_dir/csv/stem.csv."""
    csv_path = resolve_result_csv_path(output_dir / f"{stem}.csv")
    return csv_path if csv_path.exists() else None


def collect_result_csv(root: Path, file_name: str) -> list[Path]:
    """Collect CSV files and de-duplicate flat/legacy copies by logical output path."""
    if not root.exists():
        return []
    candidates = sorted(root.rglob(file_name))
    by_flat: dict[Path, Path] = {}
    for candidate in candidates:
        if candidate.suffix.lower() != ".csv":
            continue
        flat_path = normalize_result_csv_path(candidate)
        previous = by_flat.get(flat_path)
        if previous is None or previous.parent.name == "csv":
            by_flat[flat_path] = candidate
    return [by_flat[key] for key in sorted(by_flat)]


def table_preview_key(csv_path: Path) -> str:
    parts = csv_path.parts
    stem = csv_path.stem
    if "DEG" in parts and stem == "summary":
        return "deg_summary"
    if "DEG" in parts and stem == "significant_genes":
        return "deg_gene_table"
    if "run_summary" in parts and stem == "summary":
        return "tf_run_summary"
    if "TF_summary" in parts and "candidates" in parts:
        return "tf_candidates"
    if "TF_summary" in parts and "intersections" in parts and stem == "intersection_summary":
        return "tf_intersection_summary"
    if "TF_summary" in parts and "intersections" in parts and stem == "summary":
        return "tf_intersection_summary"
    if "TF_summary" in parts and "method_final" in parts and stem == "method_final_summary":
        return "tf_method_final_summary"
    if "TF_summary" in parts and "method_final" in parts and stem == "summary":
        return "tf_method_final_summary"
    if "TF_summary" in parts and "method_final" in parts:
        return "tf_method_final"
    if "intersect" in parts and stem == "summary":
        return "intersect_summary"
    if "intersect" in parts and stem == "gene_list":
        return "intersect_gene_list"
    if "intersect" in parts and stem == "deg_results":
        return "intersect_deg_results"
    if "GSEA_summary" in parts and stem == "summary":
        return "gsea_summary"
    if "GSEA" in parts and stem == "gsea_result":
        return "gsea_result"
    return "deg_gene_table"


def preview_output_name(csv_path: Path) -> str:
    try:
        relative = normalize_result_csv_path(csv_path).resolve().relative_to(RESULT_ROOT.resolve())
    except ValueError:
        relative = normalize_result_csv_path(csv_path)
    return f"{sanitize_name(relative.with_suffix('').as_posix())}.tex"


def make_preview_table(csv_path: Path, n_rows: int = TABLE_PREVIEW_ROWS) -> Path | None:
    key = table_preview_key(csv_path)
    columns = TABLE_COLUMN_PRESETS.get(key, [])
    return write_generated_latex_table(
        csv_path=resolve_result_csv_path(csv_path),
        output_name=preview_output_name(csv_path),
        columns=columns,
        n_rows=n_rows,
    )



def prefer_png_pdf(png_path: Path | None = None, pdf_path: Path | None = None) -> Path | None:
    """Prefer PNG for Beamer insertion; keep PDF as a fallback if PNG is absent."""
    if png_path is not None and png_path.exists():
        return png_path
    if pdf_path is not None and pdf_path.exists():
        return pdf_path
    return None


def first_plot_level(plot_subdir: str, fig: Path) -> str:
    """Return the analysis/configuration directory from a plot path.

    Expected layout:
      results/ngs/GSE114012/plots/<plot_subdir>/<analysis>/png/<figure>.png
    """
    parts = fig.relative_to(PLOT_ROOT / plot_subdir).parts
    return parts[0] if parts else fig.parent.name


def write_text(path: Path, lines: list[str]) -> None:
    text = "\n".join(lines) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_text(encoding="utf-8") == text:
        return
    path.write_text(text, encoding="utf-8")


def shorten_text(text: object, max_chars: int = 90) -> str:
    text = str(text)
    text = re.sub(r"\s+", " ", text).strip()
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 3] + "..."


SCI_NUMBER_PATTERN = re.compile(r"(?<![A-Za-z0-9_.+-])([+-]?(?:\d+(?:\.\d*)?|\.\d+)[eE][+-]?\d+)(?![A-Za-z0-9_.+-])")


def decimal_text_without_scientific(value: str) -> str:
    """Replace scientific-notation numbers in text with plain decimal strings."""

    def convert(match: re.Match[str]) -> str:
        token = match.group(1)
        try:
            decimal_value = Decimal(token)
        except InvalidOperation:
            return token
        plain = format(decimal_value, "f")
        if "." in plain:
            plain = plain.rstrip("0").rstrip(".")
        if plain == "-0":
            plain = "0"
        return plain

    return SCI_NUMBER_PATTERN.sub(convert, value)


def format_table_cell(value: object, max_chars: int = 90) -> str:
    if value is None:
        return "--"
    value = str(value)
    if value == "" or value.lower() == "nan":
        return "--"
    # 表格预览按字符处理；若CSV中已有科学计数法写法，仅在Beamer展示时
    # 转成普通十进制字符，避免pvalue/qvalue等字段以e-notation显示。
    value = re.sub(r"[\r\n\t]+", " ", value).strip()
    return decimal_text_without_scientific(value)


def read_csv_preview(csv_path: Path, columns: list[str], n_rows: int = 21) -> tuple[list[str], list[list[str]]]:
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        available_columns = [column for column in columns if column in (reader.fieldnames or [])]
        rows: list[list[str]] = []
        for row in reader:
            rows.append([format_table_cell(row.get(column, ""), max_chars=95) for column in available_columns])
            if len(rows) >= n_rows:
                break
    return available_columns, rows


def write_generated_latex_table(csv_path: Path, output_name: str, columns: list[str], n_rows: int = 21) -> Path | None:
    csv_path = resolve_result_csv_path(csv_path)
    if not csv_path.exists():
        return None
    headers, rows = read_csv_preview(csv_path, columns=columns, n_rows=n_rows)
    if not headers:
        return None

    tex_path = GENERATED_TABLE_ROOT / output_name
    tex_path.parent.mkdir(parents=True, exist_ok=True)
    col_spec = "|" + "|".join(["c"] * len(headers)) + "|"
    lines = [
        r"\begingroup",
        r"\fontsize{5.85pt}{7.65pt}\selectfont",
        r"\setlength{\tabcolsep}{3.35pt}",
        r"\renewcommand{\arraystretch}{1.82}",
        r"\setlength{\arrayrulewidth}{0.42pt}",
        rf"\begin{{tabular}}{{@{{}}{col_spec}@{{}}}}",
        r"\hline",
        " & ".join(r"\textbf{\textcolor{black}{" + tex_escape(header) + "}}" for header in headers) + r" \\",
        r"\hline",
    ]
    for row in rows:
        lines.append(" & ".join(tex_escape(value) for value in row) + r" \\")
        lines.append(r"\hline")
    lines.extend([
        r"\end{tabular}",
        r"\endgroup",
    ])
    write_text(tex_path, lines)
    GENERATED_PREVIEW_TABLES.add(tex_path.resolve())
    return tex_path


def make_preview_tables_parallel(csv_paths: list[Path]) -> dict[Path, Path | None]:
    """Generate many Beamer table previews in parallel."""
    unique_paths = sorted({resolve_result_csv_path(path) for path in csv_paths})
    if not unique_paths:
        return {}
    workers = min(TABLE_PREVIEW_WORKERS, len(unique_paths))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        previews = list(pool.map(make_preview_table, unique_paths))
    return dict(zip(unique_paths, previews))


def section_file(name: str) -> Path:
    return SECTION_ROOT / name


def section_cover(title: str, subtitle: str, body: str) -> list[str]:
    return [
        r"\SectionCover",
        f"  {{{tex_escape(title)}}}",
        f"  {{{tex_escape(subtitle)}}}",
        "  {",
        f"    {tex_escape(body)}",
        "  }",
        "",
    ]


def text_frame(title: str, items: list[str]) -> list[str]:
    lines = [
        f"\\TextFrame{{{tex_escape(title)}}}",
        "{",
        r"  \begin{itemize}",
        r"    \setlength{\itemsep}{0.36em}",
        r"    \setlength{\topsep}{0pt}",
        r"    \setlength{\partopsep}{0pt}",
        r"    \setlength{\parsep}{0pt}",
    ]
    lines.extend([f"    \\item {tex_escape(item)}" for item in items])
    lines.extend([
        r"  \end{itemize}",
        "}",
        "",
    ])
    return lines


def figure_frame(title: str, figure_path: Path, items: list[str], wide: bool = False) -> list[str]:
    macro = r"\ResultWideFigureFrame" if wide else r"\ResultFigureFrame"
    lines = [
        macro,
        f"  {{{tex_escape(title)}}}",
        f"  {{{rel(figure_path)}}}",
        "  {",
        r"    \begin{itemize}",
        r"      \setlength{\itemsep}{0.28em}",
        r"      \setlength{\topsep}{0pt}",
        r"      \setlength{\partopsep}{0pt}",
        r"      \setlength{\parsep}{0pt}",
    ]
    lines.extend([f"      \\item {tex_escape(item)}" for item in items])
    lines.extend([
        r"    \end{itemize}",
        "  }",
        "",
    ])
    return lines


def two_table_frame(
    title: str,
    left_tex: Path,
    left_source_csv: Path,
    left_items: list[str],
    right_tex: Path,
    right_source_csv: Path,
    right_items: list[str],
) -> list[str]:
    lines = [
        f"\\ResultTwoTableFrame{{{tex_escape(title)}}}",
        f"  {{{rel(left_source_csv)}}}",
        f"  {{{rel(left_tex)}}}",
        "  {",
        r"    \begin{itemize}",
        r"      \setlength{\itemsep}{0.24em}",
        r"      \setlength{\topsep}{0pt}",
        r"      \setlength{\partopsep}{0pt}",
        r"      \setlength{\parsep}{0pt}",
    ]
    lines.extend([f"      \\item {tex_escape(item)}" for item in left_items])
    lines.extend([
        r"    \end{itemize}",
        "  }",
        f"  {{{rel(right_source_csv)}}}",
        f"  {{{rel(right_tex)}}}",
        "  {",
        r"    \begin{itemize}",
        r"      \setlength{\itemsep}{0.24em}",
        r"      \setlength{\topsep}{0pt}",
        r"      \setlength{\partopsep}{0pt}",
        r"      \setlength{\parsep}{0pt}",
    ])
    lines.extend([f"      \\item {tex_escape(item)}" for item in right_items])
    lines.extend([
        r"    \end{itemize}",
        "  }",
        "",
    ])
    return lines


def stacked_two_table_frame(
    title: str,
    top_tex: Path,
    top_source_csv: Path,
    top_items: list[str],
    bottom_tex: Path,
    bottom_source_csv: Path,
    bottom_items: list[str],
) -> list[str]:
    lines = [
        f"\\ResultStackedTwoTableFrame{{{tex_escape(title)}}}",
        f"  {{{rel(top_source_csv)}}}",
        f"  {{{rel(top_tex)}}}",
        "  {",
        r"    \begin{itemize}",
        r"      \setlength{\itemsep}{0.18em}",
        r"      \setlength{\topsep}{0pt}",
        r"      \setlength{\partopsep}{0pt}",
        r"      \setlength{\parsep}{0pt}",
    ]
    lines.extend([f"      \\item {tex_escape(item)}" for item in top_items])
    lines.extend([
        r"    \end{itemize}",
        "  }",
        f"  {{{rel(bottom_source_csv)}}}",
        f"  {{{rel(bottom_tex)}}}",
        "  {",
        r"    \begin{itemize}",
        r"      \setlength{\itemsep}{0.18em}",
        r"      \setlength{\topsep}{0pt}",
        r"      \setlength{\partopsep}{0pt}",
        r"      \setlength{\parsep}{0pt}",
    ])
    lines.extend([f"      \\item {tex_escape(item)}" for item in bottom_items])
    lines.extend([
        r"    \end{itemize}",
        "  }",
        "",
    ])
    return lines


def table_frame(title: str, tex_path: Path, items: list[str], source_csv: Path | None = None) -> list[str]:
    source_csv = source_csv or tex_path.with_suffix(".csv")
    lines = [
        f"\\ResultTableFrame{{{tex_escape(title)}}}",
        f"  {{{rel(source_csv)}}}",
        f"  {{{rel(tex_path)}}}",
        "  {",
        r"    \begin{itemize}",
        r"      \setlength{\itemsep}{0.28em}",
        r"      \setlength{\topsep}{0pt}",
        r"      \setlength{\partopsep}{0pt}",
        r"      \setlength{\parsep}{0pt}",
    ]
    lines.extend([f"      \\item {tex_escape(item)}" for item in items])
    lines.extend([
        r"    \end{itemize}",
        "  }",
        "",
    ])
    return lines


def image_result_title(script_no: str, analysis: str, result_label: str) -> str:
    """Create a concise slide title that states the analysis content."""
    section_title = SCRIPT_TITLES.get(script_no, (script_no, ""))[0]
    return f"{script_no}. {section_title} | {analysis} | {result_label}"


def deg_table_title(analysis: str, stem: str) -> str:
    labels = {
        "summary": "DEG统计摘要",
        "significant_genes": "显著差异基因表",
    }
    return f"01. 差异表达分析 | {analysis} | {labels.get(stem, stem)}"


def intersect_table_title(scheme: str, rel_parts: tuple[str, ...]) -> str:
    file_name = rel_parts[-1]
    if file_name == "summary.csv":
        return f"02. 显著基因交集 | {scheme} | 交集统计摘要"
    if file_name == "gene_list.csv":
        return f"02. 显著基因交集 | {scheme} | 交集基因注释列表"
    analysis = rel_parts[0] if len(rel_parts) > 1 else scheme
    return f"02. 显著基因交集 | {scheme} | {analysis}中的交集基因DEG结果"


def gsea_table_title(analysis: str, geneset: str) -> str:
    return f"06. GSEA富集结果 | {analysis} | {gsea_display_label(geneset)}"


def gsea_pair_title(analysis: str, geneset: str, result_type: str) -> str:
    return f"06/07. GSEA图表配对 | {analysis} | {gsea_display_label(geneset)} | {result_type}"


def tf_table_title(input_type: str, analysis: str, parts: tuple[str, ...]) -> str:
    if parts == ("method_final", "summary", "method_final_summary.csv"):
        return f"09. TF整合 | {input_type}/{analysis} | 六种方法final结果摘要"
    if parts == ("intersections", "summary", "intersection_summary.csv"):
        return f"09. TF整合 | {input_type}/{analysis} | TF交集方案总览"
    if len(parts) >= 4 and parts[0] == "intersections":
        scheme = parts[1]
        if parts[2] == "summary":
            return f"09. TF整合 | {input_type}/{analysis} | {tf_intersection_label(scheme)} | 交集摘要"
        if parts[2] == "candidates":
            return f"09. TF整合 | {input_type}/{analysis} | {tf_intersection_label(scheme)} | Top10候选TF"
    if len(parts) >= 3 and parts[0] == "method_final":
        method = parts[1]
        return f"09. TF整合 | {input_type}/{analysis} | {tf_method_label(method)} final TF排序"
    return f"09. TF整合 | {input_type}/{analysis} | {' / '.join(parts)}"


def compact_table_frame(title: str, rows: list[tuple[str, str, str]]) -> list[str]:
    lines = [
        f"\\begin{{frame}}{{{tex_escape(title)}}}",
        r"  \begin{adjustbox}{max width=0.96\textwidth,max totalheight=0.78\textheight}",
        r"  \begin{tabular}{@{}lll@{}}",
        r"    \toprule",
        r"    编号 & 脚本 & 当前展示重点 \\",
        r"    \midrule",
    ]
    for code, script, focus in rows:
        lines.append(f"    {tex_escape(code)} & {tex_escape(script)} & {tex_escape(focus)} \\\\")
    lines.extend([
        r"    \bottomrule",
        r"  \end{tabular}",
        r"  \end{adjustbox}",
        r"\end{frame}",
        "",
    ])
    return lines


def build_project_design() -> str:
    section_name = "project/project_design.tex"
    lines: list[str] = []
    lines += section_cover(
        "课题设计：结直肠癌休眠细胞苏醒与复发转移",
        "从远期复发临床问题出发，解析休眠癌细胞苏醒的转录调控与可干预节点",
        "结直肠癌远处复发可能并不只是癌细胞是否播散的问题，更关键的是长期休眠的播散性肿瘤细胞是否被重新激活。当前分析围绕休眠样 LRC 与循环样 BULK 的转录差异，逐步提炼稳定候选基因、通路和转录因子。",
    )
    lines += text_frame(
        "研究背景与核心问题",
        [
            "结直肠癌是全球最常见的恶性肿瘤之一。临床上一个长期存在但尚未完全解决的问题是，许多患者在接受根治性手术后数年甚至十余年仍会发生远处转移复发。",
            "越来越多的证据表明，部分结直肠癌细胞可能在疾病早期便已完成播散，并定植于肝脏、肺脏或骨髓等远处器官。",
            "然而，这些播散性肿瘤细胞（Disseminated Tumor Cells, DTCs）并不会立即形成转移灶，而是长期处于一种低增殖甚至静止的休眠状态。",
            "这一现象提示，对于部分患者而言，“转移的发生”可能早已在疾病早期被决定，而真正决定患者长期预后的关键因素并非癌细胞是否播散，而是这些休眠癌细胞是否被重新激活和苏醒。",
            "当休眠细胞受到炎症反应、组织重塑、免疫微环境改变或代谢重编程等因素刺激后，可重新进入细胞周期并形成临床可检测的转移灶，最终导致疾病复发。",
        ],
    )
    lines += text_frame(
        "研究目的与转化价值",
        [
            "因此，与传统研究聚焦于抑制肿瘤转移不同，维持播散癌细胞的长期休眠状态可能是一种更具临床转化价值的治疗策略。",
            "本研究拟系统解析结直肠癌休眠细胞向苏醒细胞转变过程中的关键转录因子调控网络、信号通路变化及免疫微环境重塑特征。",
            "并进一步通过药物重定位策略筛选能够维持肿瘤休眠状态的候选药物，为预防结直肠癌复发和转移提供新的理论依据和治疗靶点。",
        ],
    )
    lines += text_frame(
        "整体研究步骤",
        [
            "步骤一：定义研究模型。找结直肠癌 dormancy / recurrence / minimal residual disease / liver metastasis relapse 数据集。优先分组：原发癌、转移癌、复发癌、治疗后残留癌细胞、休眠样癌细胞。",
            "步骤二：建立休眠评分。用 dormancy markers 计算休眠 score。",
            "步骤三：识别“苏醒癌细胞”。单细胞中按 dormancy score 和 proliferation score 分群：休眠型、过渡型、苏醒型、增殖型。核心比较：休眠型 vs 苏醒型。",
            "步骤四：寻找导致苏醒的转录因子。做 SCENIC / DoRothEA / pySCENIC / ChIP-X enrichment。",
            "步骤五：构建苏醒轨迹。用 Monocle3 / Slingshot / CytoTRACE。构建：休眠 → 过渡 → 苏醒 → 增殖。沿 pseudotime 找动态 TF 和动态通路。",
            "步骤六：解析苏醒相关通路。做 GSVA / GSEA。",
            "步骤七：分析免疫微环境。看苏醒型癌细胞是否伴随 M2 macrophage、TREM2/APOE macrophage、CAF、Treg、exhausted CD8 增加。用 CellChat / NicheNet 分析癌细胞与免疫/基质细胞互作。",
            "步骤八：提炼关键调控轴。例如：CAF-derived TGFβ/IL6 → STAT3/YAP → 癌细胞苏醒。Macrophage TNF/IL1β → NFκB → 癌细胞苏醒。WNT niche → β-catenin/TCF → 癌细胞苏醒。",
            "步骤九：Drug repurposing。将“苏醒型上调基因”和“关键 TF 靶基因”输入 CMap / LINCS / DGIdb / Enrichr Drug Signatures。筛选能逆转苏醒 signature、维持休眠状态的候选药物。Alphafold+分子对接证明直接结合。",
            "步骤十：形成最终假说。结直肠癌复发不是单纯因为癌细胞增殖，而是微环境信号诱导休眠癌细胞苏醒。阻断关键 TF/通路，可让残留癌细胞长期维持休眠，从而减少复发和转移。",
        ],
    )
    write_text(section_file(section_name), lines)
    return section_name


def build_00() -> str:
    section_name = f"{DATASET_ID}/00_sample_clustering_heatmap.tex"
    lines = section_cover(
        "00. 样本结构质控：LRC/BULK 是否具备可比较性",
        "用 TPM 相关性和层级聚类检查休眠样模型的样本基础",
        "在进入 DEG、GSEA 和 TF 推断之前，必须先确认样本在表达层面没有明显混乱或离群。00 号结果用于回答：LRC 与 BULK 样本是否形成可解释的表达结构，后续比较是否具有可信起点。",
    )
    lines += text_frame(
        "为什么先做样本聚类",
        [
            "输入使用 SE 对象中的 TPM，而不是原始 count；这样样本间相关性更适合表达模式质控。",
            "样本标签使用 clinical table 的 Title 列，便于直接识别细胞系、重复和 LRC/BULK 状态。",
            "聚类只由表达相关性决定，不强制把同一细胞系或同一分组放在一起；如果图中出现离群或异常聚类，后续 DEG 和富集结果需要谨慎解释。",
            "该步骤在证据链中的作用是确认模型结构：我们后续讨论的休眠样转录程序，必须建立在样本表达结构可解释的基础上。",
        ],
    )
    plot_dir = PLOT_ROOT / "sample_clustering_heatmap"
    for fig in sorted(plot_dir.rglob("heatmap.png")):
        parts = fig.relative_to(plot_dir).parts
        label = " / ".join(parts[:-2]) if len(parts) > 2 else fig.parent.name
        lines += figure_frame(
            image_result_title("00", label, "样本相关性层级聚类"),
            fig,
            [
                "热图主体为样本间 TPM 相关性矩阵；颜色越接近高相关端，说明两个样本整体表达谱越相似。",
                "层级聚类树反映表达模式驱动的样本相似性，而不是预设分组；因此可用于识别潜在离群样本、批次结构或细胞系特异性。",
                "若 LRC 或 BULK 样本在该图中呈现相对稳定结构，后续 LRC vs BULK 的差异表达、GSEA 和 TF 推断才更容易被解释为生物状态差异，而非样本质量问题。",
            ],
        )
    write_text(section_file(section_name), lines)
    return section_name


def build_01() -> str:
    section_name = f"{DATASET_ID}/01_limma_differential_expression.tex"
    lines = section_cover(
        "01. 差异表达证据链：每个 LRC/BULK 设计的统计、基因、火山图与表达热图",
        "同一分析设计内部按 summary → significant genes → volcano → Top DEG heatmap 顺序展示",
        "01 号差异分析用于定义 LRC 休眠样状态相对 BULK 循环样状态的转录改变。这里按分析设计逐套展示，目的是让每个模型都能同时看到统计规模、显著基因明细、效应量分布和样本层面表达分离。",
    )
    lines += text_frame(
        "差异分析结果的技术含义与阅读顺序",
        [
            "summary 先回答“这个比较强不强”：Total_Significant_Genes 必须等于 Up 与 Down 的总和；Up/Down 方向定义为 LRC 相对 BULK。",
            "significant_genes 展示通过 01 号脚本阈值的基因。logFC 表示 LRC 方向效应量，t/P.Value/adj.P.Val 表示统计证据，Symbol/Ensembl/Entrez 用于后续富集和验证。",
            "传统火山图把显著性与 logFC 同时展示，适合快速判断该分析设计中显著基因是否呈现方向偏倚、强效候选是否集中，以及是否存在明显上下调不平衡。",
            "Top DEG 表达热图回到样本层面检查这些候选基因是否真正把 LRC 与 BULK 分开；如果统计显著但热图不能分离样本，后续机制解读应更谨慎。",
            "全量差异基因排序表仍作为 GSEA 输入保存在结果目录，但不再进入 Beamer 展示，避免报告被中间排序表稀释。",
        ],
    )
    analysis_dirs = sorted(
        [analysis_dir for analysis_dir in dirs(TABLE_ROOT) if (analysis_dir / "DEG").exists()],
        key=lambda path: sort_analysis_name(path.name),
    )
    for analysis_dir in analysis_dirs:
        deg_dir = analysis_dir / "DEG"
        analysis = analysis_dir.name
        lines += text_frame(
            f"01. {analysis}：LRC vs BULK 差异分析展示顺序",
            [
                f"本组展示 {analysis} 分析设计。先看 summary 判断显著 DEG 总量和上下调方向；再看 significant_genes 确认候选基因及统计量。",
                "随后用传统火山图检查 logFC 与显著性分布，最后用 Top DEG 热图验证这些差异基因能否在样本层面形成 LRC/BULK 分离。",
                "这种顺序把“统计是否成立”“候选是谁”“方向是否清楚”“样本是否支持”四个问题放在同一分析设计内连续回答。",
            ],
        )

        summary_csv = find_result_csv(deg_dir, "summary")
        if summary_csv is not None:
            tex = make_preview_table(summary_csv)
            if tex is not None:
                lines += table_frame(
                    deg_table_title(analysis, "summary"),
                    tex,
                    [
                        "该 summary 是当前分析设计的质量控制入口。重点核对 Up、Down、Total_Significant_Genes 和阈值列；Total_Significant_Genes 应等于 Up+Down。",
                        "若显著基因总量过少，说明该模型的 LRC/BULK 差异较弱；若上下调极度不平衡，后续解释要结合火山图和热图确认是否由少数强效基因驱动。",
                    ],
                    source_csv=summary_csv,
                )

        significant_csv = find_result_csv(deg_dir, "significant_genes")
        if significant_csv is not None:
            tex = make_preview_table(significant_csv)
            if tex is not None:
                lines += table_frame(
                    deg_table_title(analysis, "significant_genes"),
                    tex,
                    [
                        "该表为通过当前阈值的显著 DEG 前 21 行预览。logFC>0 表示 LRC 方向上调，logFC<0 表示 BULK 方向更高。",
                        "t、P.Value 和 adj.P.Val 用于评估统计强度；Symbol/Ensembl/Entrez 是后续交集、GSEA leading-edge 回查和 TF 富集输入时最常用的标识。",
                    ],
                    source_csv=significant_csv,
                )

        volcano_fig = PLOT_ROOT / "volcano" / analysis / "png" / "volcano_plot.png"
        if volcano_fig.exists():
            lines += figure_frame(
                image_result_title("01", analysis, "传统火山图"),
                volcano_fig,
                [
                    "该图横坐标为 logFC，纵坐标为显著性指标；红色 Sig_Up 与蓝色 Sig_Down 分别对应 LRC 方向上调和下调显著基因。",
                    "阅读时先看红蓝点数量和左右分布，再看被标注的候选基因是否位于高显著性、高效应量区域。",
                    "如果某些细胞系火山图中信号很弱，说明该模型对休眠样差异程序的贡献有限；如果多个模型均出现相似强信号，则更支持跨模型稳定候选。",
                ],
            )

        heatmap_fig = PLOT_ROOT / "gene_heatmap" / analysis / "png" / "gene_heatmap.png"
        if heatmap_fig.exists():
            lines += figure_frame(
                image_result_title("01", analysis, "Top DEG表达热图"),
                heatmap_fig,
                [
                    "该图选取当前分析设计中排序靠前的 Top DEG，基于表达矩阵进行行标准化展示；左侧方向条保留 Up/Down 信息，顶部条带显示 Group 与 Cell_Line。",
                    "阅读重点是 LRC 样本是否在这些基因上形成相对一致的表达模式，并与 BULK 样本分离；这比只看 p 值更能体现候选基因的样本层面可解释性。",
                    "如果热图中某些基因同时具有清楚分组分离、强 logFC，并出现在 02 号交集或 09 号 TF 证据链中，应优先作为后续验证对象。",
                ],
                wide=True,
            )
    write_text(section_file(section_name), lines)
    return section_name


def build_02() -> str:
    section_name = f"{DATASET_ID}/02_intersect_significant_genes.tex"
    lines = section_cover(
        "02. 显著 DEG 交集：从多模型中提炼稳定候选基因",
        "每个交集方案按 summary → gene list → 成员 DEG 结果 → 多组火山图顺序展示",
        "02 号脚本把多个 LRC/BULK 分析设计中的 significant_genes 取交集，用于降低单一细胞系背景噪音。当前展示逻辑强调：先看交集规模，再看候选基因注释，再回查每个成员 DEG 的统计证据，最后用多组火山图观察组合模型方向。",
    )
    lines += text_frame(
        "交集结果的技术含义与阅读顺序",
        [
            "summary 回答交集是否足够严格：Total_Intersected_Genes 是交集规模；Common_Up/Common_Down 表示方向一致候选；Mixed_Direction 提示不同模型方向不完全一致。",
            "gene_list 是只含基因注释的清单，用于后续 TF 富集、人工筛选或实验验证设计；它不展示 p 值和 logFC，因为这些统计证据需要回到成员 DEG 结果中查看。",
            "每个成员分析的 deg_results 展示交集基因在该 DEG 结果中的统计量。若同一个 Symbol 在多个成员中方向一致且 p 值稳定，可信度高于只在单模型出现的基因。",
            "多组火山图用于把参与交集的分析设计放在同一图中对比，判断该交集方案的显著基因方向是否具有一致趋势。",
        ],
    )
    for scheme_dir in sorted(dirs(INTERSECT_ROOT), key=lambda path: sort_analysis_name(path.name)):
        scheme = scheme_dir.name
        lines += text_frame(
            f"02. {scheme}：交集方案展示顺序",
            [
                f"本组展示 {scheme} 交集方案。先看 summary 判断交集规模与方向一致性；再看 gene_list 确认候选基因注释。",
                "随后逐个展示参与交集的 DEG 结果，回查每个交集基因在原始分析中的 logFC 与 p 值；最后展示对应多组火山图，观察不同成员分析的整体方向分布。",
            ],
        )

        summary_csv = find_result_csv(scheme_dir, "summary")
        if summary_csv is not None:
            tex = make_preview_table(summary_csv)
            if tex is not None:
                lines += table_frame(
                    intersect_table_title(scheme, ("summary.csv",)),
                    tex,
                    [
                        "该 summary 用于判断当前交集方案是否足够严格且仍保留可解释候选。重点看 Total_Intersected_Genes、Common_Up、Common_Down 与 Mixed_Direction。",
                        "如果 Mixed_Direction 较多，说明候选基因在不同模型间方向不稳定，后续不应仅凭是否进入交集来判断其生物学可靠性。",
                    ],
                    source_csv=summary_csv,
                )

        gene_list_csv = find_result_csv(scheme_dir, "gene_list")
        if gene_list_csv is not None:
            tex = make_preview_table(gene_list_csv)
            if tex is not None:
                lines += table_frame(
                    intersect_table_title(scheme, ("gene_list.csv",)),
                    tex,
                    [
                        "该 gene_list 是交集基因注释列表，只保留 Symbol/Ensembl/Entrez 等可追踪标识，适合作为后续 TF 富集、候选验证和文献检索的输入。",
                        "这里不展示统计量，是为了避免把不同 DEG 结果的 p 值混在一个清单里；统计证据需要在后续成员 DEG 结果页逐个回查。",
                    ],
                    source_csv=gene_list_csv,
                )

        member_deg_files = [
            csv_path
            for csv_path in collect_result_csv(scheme_dir, "deg_results.csv")
            if normalize_result_csv_path(csv_path).name == "deg_results.csv"
        ]
        member_deg_files = sorted(member_deg_files, key=lambda path: sort_analysis_name(normalize_result_csv_path(path).parent.name))
        for csv_path in member_deg_files:
            flat_csv = normalize_result_csv_path(csv_path)
            rel_parts = flat_csv.relative_to(scheme_dir).parts
            member = rel_parts[0] if len(rel_parts) > 1 else scheme
            tex = make_preview_table(csv_path)
            if tex is None:
                continue
            lines += table_frame(
                intersect_table_title(scheme, rel_parts),
                tex,
                [
                    f"该表展示 {scheme} 交集基因在 {member} 原始 DEG 结果中的统计证据。重点看 Symbol、logFC、t、P.Value 与 adj.P.Val。",
                    "如果同一基因在多个成员分析中 logFC 方向一致，并且 p 值稳定，则更适合作为跨模型稳定的休眠样候选基因。",
                ],
                source_csv=csv_path,
            )

        multiple_volcano_fig = PLOT_ROOT / "multiple_volcano" / scheme / "png" / "multiple_volcano_plot.png"
        if multiple_volcano_fig.exists():
            lines += figure_frame(
                image_result_title("02", scheme, "多组火山图"),
                multiple_volcano_fig,
                [
                    "该图把当前交集方案中的多个 LRC/BULK 分析设计并列展示，红色 Sig_Up 与蓝色 Sig_Down 对应各模型显著上/下调基因。",
                    "阅读时先看不同组之间上下调基因分布是否相似，再看中心组名区域两侧是否存在稳定的强效差异基因。",
                    "若多组火山图中多个成员均显示相似方向与强度，说明该交集方案更可能捕捉到跨模型休眠样转录程序，而不仅是单细胞系噪音。",
                ],
                wide=True,
            )
    write_text(section_file(section_name), lines)
    return section_name


def build_06() -> list[str]:
    section_names = [f"{DATASET_ID}/06_gsea_analysis.tex"]
    lines = section_cover(
        "06. GSEA 运行概览：把 DEG 排序转化为通路方向",
        "基于全量 ranked genes 的 MSigDB 富集，不只依赖显著基因列表",
        "06 号脚本把每套 LRC/BULK DEG 的全量排序基因输入 clusterProfiler::GSEA，并按指定 MSigDB 类别批量输出结果。GSEA 用来回答：休眠样 LRC 的表达改变是否集中在特定通路、免疫/炎症程序、转录因子靶集或肿瘤相关 signature 上。",
    )
    lines += text_frame(
        "GSEA 结果如何服务于休眠复苏假说",
        [
            "GSEA 使用 all_genes 的全量排序，当前 rank metric 为 t statistic；这比只看 significant_genes 更适合捕捉整体通路偏移。",
            "06 号脚本负责运算和保存 GSEA 结果，07 号脚本负责绘图；Beamer 中不再展示全部 analysis × gene set 类别的运行清单，避免用流程性 summary 稀释结果展示。",
            "当前展示只保留已经生成 dotplot 与 gsea_result.csv 的组合，并按同一 GSEA 设计先图后表展示，便于把视觉结论与统计字段逐一对应。",
            "正 NES 通常提示 LRC 方向富集，负 NES 通常提示 BULK 方向富集；这为“休眠维持”或“苏醒启动”相关通路筛选提供方向性。",
            "阅读时优先关注与炎症反应、组织重塑、免疫调节、细胞周期、代谢重编程、TFT 靶集或肿瘤相关 signature 相关的通路。",
        ],
    )
    write_text(section_file(section_names[0]), lines)

    return section_names


def build_07() -> list[str]:
    section_names = [f"{DATASET_ID}/07_gsea_plotting.tex"]
    lines = section_cover(
        "07. GSEA 图表配对解读：先看通路图，再看统计表",
        "每个 analysis/geneset 的 dotplot 与同一结果表连续展示",
        "本章节把每个 GSEA 设计的图片和表格放在一起阅读。dotplot 用于快速识别 top 通路，紧随其后的表格用于核对 NES、pvalue、p.adjust 和 qvalue 等核心统计字段。",
    )
    lines += text_frame(
        "GSEA 图表配对的阅读顺序",
        [
            "先看 dotplot：确认哪些通路进入 top10，以及通路名称是否与休眠、炎症、ECM remodeling、代谢重编程、免疫微环境或 TF 靶集相关。",
            "再看结果表：本报告只展示 ID、NES、pvalue、p.adjust、qvalue 这些最关键字段；NES 判断富集方向，p.adjust/qvalue 判断多重校正后的统计证据。",
            "07 号 dotplot 的显著性筛选与 06 号 GSEA 运算保持一致：当前统一使用 p.adjust < 0.05。",
            "单通路 GSEA+热图数量较多，暂不纳入逐页汇报；需要深入某一通路时可回到结果目录查看对应单通路图。",
        ],
    )
    write_text(section_file(section_names[0]), lines)

    by_analysis: dict[str, list[Path]] = {}
    for csv_path in collect_result_csv(TABLE_ROOT, "gsea_result.csv"):
        flat_csv = normalize_result_csv_path(csv_path)
        parts = flat_csv.relative_to(TABLE_ROOT).parts
        if len(parts) < 4 or parts[1] != "GSEA":
            continue
        analysis = parts[0]
        by_analysis.setdefault(analysis, []).append(csv_path)

    for analysis, csv_files in sorted(by_analysis.items(), key=lambda item: sort_analysis_name(item[0])):
        name = f"{DATASET_ID}/07_gsea_paired_results_{sanitize_name(analysis)}.tex"
        section_names.append(name)
        lines = section_cover(
            f"07. GSEA paired results：{analysis}",
            f"{analysis} 的 GSEA 图表配对结果",
            "每个 MSigDB 类别先展示 dotplot 图，再展示同一 analysis/geneset 对应的 GSEA 结果表。图中 top 通路可直接与下一页 NES、pvalue、p.adjust 和 qvalue 对应查看。",
        )
        csv_files = sorted(csv_files, key=sort_gsea_csv)
        preview_map = make_preview_tables_parallel(csv_files)
        for csv_path in csv_files:
            flat_csv = normalize_result_csv_path(csv_path)
            geneset = flat_csv.relative_to(TABLE_ROOT).parts[2]
            fig = PLOT_ROOT / "GSEA" / analysis / geneset / "png" / "dotplot.png"
            if fig.exists():
                lines += figure_frame(
                    gsea_pair_title(analysis, geneset, "Dotplot"),
                    fig,
                    [
                        "该 dotplot 展示当前 analysis × MSigDB 类别中通过 p.adjust < 0.05 后的 top10 通路，用于快速定位最强富集主题。",
                        "气泡大小和颜色用于同时表达富集规模和统计显著性；通路名称若与休眠、炎症、免疫、细胞周期、ECM 或 TF 靶集相关，应优先进入后续解释。",
                        "下一页是同一 GSEA 设计的结果表，请用 NES、p.adjust 与 qvalue 验证图中通路是否具有可靠统计证据和明确富集方向。",
                    ],
                    wide=True,
                )

            compact_tex = preview_map.get(resolve_result_csv_path(csv_path))
            if compact_tex is None:
                continue
            lines += table_frame(
                gsea_pair_title(analysis, geneset, "GSEA结果表"),
                compact_tex,
                [
                    "该表是同一 dotplot 对应的 GSEA 统计结果预览；为保证汇报版式，只展示 ID、NES、pvalue、p.adjust、qvalue 等核心判断字段。",
                    "NES 是标准化富集分数：正值通常解释为 LRC 方向富集，负值通常解释为 BULK 方向富集。",
                    "p.adjust 是本报告与 06 号运算统一采用的显著性判断字段，当前阈值为 p.adjust < 0.05；qvalue 可作为额外稳健性参考。",
                    "如果某个通路在 dotplot 中显著且 NES 方向明确，后续可回到完整 CSV 查看被隐藏的 leading_edge、rank、Description 等详细字段。",
                ],
                source_csv=csv_path,
            )
        write_text(section_file(name), lines)
    return section_names


def tf_input_key(path: Path, root: Path) -> tuple[str, str]:
    parts = path.relative_to(root).parts
    if len(parts) >= 2:
        return parts[0], parts[1]
    return "unknown", "unknown"


def build_08() -> list[str]:
    section_names = [f"{DATASET_ID}/08_tf_enrichment_analysis.tex"]
    lines = section_cover(
        "08. TF 富集/活性推断：寻找休眠样转录程序的上游调控因子",
        "六种方法分别从基因列表、调控网络、数据库证据和 signed regulon 角度推断 TF",
        "08 号脚本产生的是原始 TF 证据；不同方法的输入、统计量和方向性并不完全一致，因此本报告不逐页展示原始表，而是把 08 号输出交给 09 号统一整理、排序和取交集。",
    )
    lines += text_frame(
        "六种 TF 方法在证据链中的定位",
        [
            "DoRothEA 与 TRRUST 更偏向 TF-target overlap，适合快速筛选与显著基因集合存在靶基因重叠的 TF。",
            "ChEA3 与 ENRICHR 提供多数据库/多 library 的排名证据，适合判断某个 TF 是否在多个外部证据源中反复出现。",
            "VIPER 与 CollecTRI 使用带方向的 regulon 或 signed network，可以提供 TF activity 倾向，更适合解释“哪个 TF 可能在 LRC 中被激活或抑制”。",
            "输入分两类：DEG significant_genes 用于单个/组合分析设计，intersect gene_list 用于跨模型稳定候选；这两类输入分别回答模型特异和跨模型稳定 TF 线索。",
            "09 号整合脚本会把这些异构结果压缩成 method_final、intersection_summary 和 Top10 candidate 表，便于最终人工筛选。",
        ],
    )
    write_text(section_file(section_names[0]), lines)
    return section_names


def build_09() -> list[str]:
    section_names = [f"{DATASET_ID}/09_integrate_tf_enrichment_results.tex"]
    lines = section_cover(
        "09. TF 候选整合：从六方法结果中提炼可验证调控因子",
        "按 DEG 与 intersect 输入方案分别展示 method_final、交集方案摘要和 Top10 TF 候选",
        "09 号结果是 TF 分析的汇报重点。它把六种方法的原始输出整理成可横向比较的 final 排名，并通过多种交集方案筛选更稳健的候选 TF，用于后续构建休眠维持或苏醒启动调控轴。",
    )
    lines += text_frame(
        "TF 整合结果的阅读顺序",
        [
            "每个输入方案对应一个目录层级：DEG/<analysis> 表示来自某套差异分析的显著基因；intersect/<scheme> 表示来自 02 号跨模型交集基因。",
            "先看 TF 交集方案摘要，判断哪些方法组合得到的候选数量更少、更严格；再看六种方法 final 结果，理解每种方法给出的排序基础。",
            "最后看 Top10 候选交集表。Consensus_Rank 越靠前、Source_Method_Count 越高、CheA3_Library_Count 越大，通常说明该 TF 更值得进入后续验证。",
            "需要注意：不同方法的统计量含义不同。VIPER/CollecTRI 的 NES/Direction 更接近活性方向；DoRothEA/TRRUST/ENRICHR/ChEA3 更偏向靶基因或数据库证据支持。",
        ],
    )
    summary_csv = find_result_csv(TF_SUMMARY_ROOT / "run_summary", "summary")
    if summary_csv is not None:
        summary_tex = make_preview_table(summary_csv)
        if summary_tex is not None:
            lines += table_frame(
                "09. TF整合运行总览 | 全部输入方案",
                summary_tex,
                [
                    "该表概括 09 号脚本整合了多少 TF 输入方案和交集方案。Input_Type 区分 DEG 与 intersect；TF_Analysis_Name 对应差异分析设计或交集方案。",
                    "Methods_Integrated 表示纳入的 TF 方法数量；Intersection_Schemes 表示后续会展示多少种多方法交集组合。",
                ],
                source_csv=summary_csv,
            )
    write_text(section_file(section_names[0]), lines)

    all_csv = collect_result_csv(TF_SUMMARY_ROOT, "*.csv")
    grouped: dict[tuple[str, str], list[Path]] = {}
    for csv_path in all_csv:
        if "/run_summary/" in csv_path.as_posix():
            continue
        rel_parts = normalize_result_csv_path(csv_path).relative_to(TF_SUMMARY_ROOT).parts
        if len(rel_parts) < 2:
            continue
        input_type, analysis = rel_parts[0], rel_parts[1]
        grouped.setdefault((input_type, analysis), []).append(csv_path)
    for (input_type, analysis), csv_files in sorted(grouped.items(), key=lambda item: sort_tf_group_key(item[0][0], item[0][1])):
        name = f"{DATASET_ID}/09_tf_summary_{sanitize_name(input_type)}__{sanitize_name(analysis)}.tex"
        section_names.append(name)
        input_label = input_type.upper() if input_type.lower() == "deg" else "intersect"
        lines = section_cover(
            f"09. TF结果整合：{input_label} / {analysis}",
            "按目录结构展示该输入方案的 TF final 结果和多方法交集候选",
            "本分文件对应 TF_summary 下的一个输入方案。先看交集方案摘要，再看六种方法各自 final 排序，最后看多方法交集 Top10 候选。这样可以同时保留单方法证据和跨方法稳定性。",
        )
        csv_files = sorted(csv_files, key=sort_tf_csv)
        preview_map = make_preview_tables_parallel(csv_files)

        method_summary: Path | None = None
        method_results: list[Path] = []
        intersection_summary: Path | None = None
        scheme_summary: dict[str, Path] = {}
        scheme_candidates: dict[str, Path] = {}

        for csv_path in csv_files:
            flat_csv = normalize_result_csv_path(csv_path)
            parts = flat_csv.relative_to(TF_SUMMARY_ROOT / input_type / analysis).parts
            if parts == ("method_final", "summary", "method_final_summary.csv"):
                method_summary = csv_path
            elif len(parts) >= 3 and parts[0] == "method_final" and parts[1] != "summary":
                method_results.append(csv_path)
            elif parts == ("intersections", "summary", "intersection_summary.csv"):
                intersection_summary = csv_path
            elif len(parts) >= 4 and parts[0] == "intersections" and parts[2] == "summary":
                scheme_summary[parts[1]] = csv_path
            elif len(parts) >= 4 and parts[0] == "intersections" and parts[2] == "candidates":
                scheme_candidates[parts[1]] = csv_path

        if intersection_summary is not None:
            tex = preview_map.get(resolve_result_csv_path(intersection_summary))
            if tex is not None:
                lines += table_frame(
                    f"09. TF交集方案摘要 | {input_label}/{analysis}",
                    tex,
                    [
                        "该页先展示当前输入方案下全部 TF 交集方案摘要，用于判断不同组合的严格程度和候选规模。",
                        "Intersection_Name 为交集方案名；Required_Methods 列出参与方法；Number_Of_Methods 为方法数量；Intersected_TF_Count 为全量交集 TF 数；Reported_Top_N 为后续候选表展示数量。",
                        "交集数量越少通常越严格；若严格交集仍保留多个可解释 TF，说明该输入方案具有较强跨方法一致性。",
                    ],
                    source_csv=intersection_summary,
                )

        if method_summary is not None:
            tex = preview_map.get(resolve_result_csv_path(method_summary))
            if tex is not None:
                lines += table_frame(
                    f"09. 六方法整合摘要 | {input_label}/{analysis}",
                    tex,
                    [
                        "该页展示六种 TF 富集/活性推断方法整理后的 final TF 数量，是后续多方法交集的输入来源。",
                        "Method 表示方法名；Final_TF_Count 表示该方法最终进入排序和交集分析的 TF 数目。",
                        "如果某方法 Final_TF_Count 很低，说明该方法在当前基因集或 ranked signature 上覆盖有限，解释交集结果时需要考虑方法覆盖度。",
                    ],
                    source_csv=method_summary,
                )

        for csv_path in sorted(method_results, key=sort_tf_csv):
            flat_csv = normalize_result_csv_path(csv_path)
            parts = flat_csv.relative_to(TF_SUMMARY_ROOT / input_type / analysis).parts
            method = parts[1]
            tex = preview_map.get(resolve_result_csv_path(csv_path))
            if tex is None:
                continue
            lines += table_frame(
                f"09. 单方法TF整合结果 | {input_label}/{analysis} | {tf_method_label(method)}",
                tex,
                [
                    f"该页展示 {tf_method_label(method)} 方法整理后的 final TF 排名。该表用于回答：在单一方法证据下，哪些 TF 最可能解释当前输入基因集。",
                    "Rank 为该方法内部排序；TF 为转录因子 symbol；Score/P_Value/Adjusted_P_Value/NES/Direction 为该方法可提供的统计量或活性方向。",
                    "CheA3_Library_Count 表示该 TF 在 ChEA3 多 library 中出现的证据数量；CheA3_Integrated_TopRank 用于补充跨数据库可靠性判断。",
                    "单方法排名不是最终结论，后续需要结合交集候选表判断其跨方法稳定性。",
                ],
                source_csv=csv_path,
            )

        ordered_schemes = sorted(
            set(scheme_candidates),
            key=lambda scheme: (TF_INTERSECTION_ORDER.get(scheme, 99), scheme),
        )
        candidate_entries: list[tuple[str, Path, Path]] = []
        for scheme in ordered_schemes:
            csv_path = scheme_candidates.get(scheme)
            if csv_path is None:
                continue
            tex = preview_map.get(resolve_result_csv_path(csv_path))
            if tex is not None:
                candidate_entries.append((scheme, tex, csv_path))

        candidate_note = [
            "Consensus_Rank 为综合排序；TF 为转录因子 symbol；Source_Methods 显示该候选来自哪些方法组合。",
            "Source_Method_Count 表示支持该 TF 的方法数；CheA3_Library_Count 表示 ChEA3 多证据库支持数，数量越高通常说明外部证据越丰富。",
            "DoRothEA/ChEA3/VIPER/ENRICHR/TRRUST/CollecTRI 列保留各方法紧凑证据。优先关注多方法支持、ChEA3 library 较多、且在活性方法中方向清楚的 TF。",
        ]
        for i in range(0, len(candidate_entries), 2):
            first = candidate_entries[i]
            second = candidate_entries[i + 1] if i + 1 < len(candidate_entries) else None
            first_title = tf_intersection_label(first[0])
            first_items = [f"上表：{first_title}。"] + candidate_note
            if second is None:
                lines += table_frame(
                    f"09. TF交集Top10候选 | {input_label}/{analysis} | {first_title}",
                    first[1],
                    first_items,
                    source_csv=first[2],
                )
                continue
            second_title = tf_intersection_label(second[0])
            second_items = [f"下表：{second_title}。"] + candidate_note
            lines += stacked_two_table_frame(
                f"09. TF交集Top10候选 | {input_label}/{analysis} | {first_title} / {second_title}",
                first[1],
                first[2],
                first_items,
                second[1],
                second[2],
                second_items,
            )
        write_text(section_file(name), lines)
    return section_names


def existing_tf_summary_dir(input_type: str, analysis: str) -> Path | None:
    """Find a TF_summary directory while tolerating DEG/deg naming variants."""
    for candidate_type in (input_type, input_type.lower(), input_type.upper()):
        candidate = TF_SUMMARY_ROOT / candidate_type / analysis
        if candidate.exists():
            return candidate
    return None


def collect_tf_summary_csv(input_type: str, analysis: str) -> list[Path]:
    """Collect generated 09号TF整合结果 for one DEG/intersect input scheme."""
    summary_dir = existing_tf_summary_dir(input_type, analysis)
    if summary_dir is None:
        return []
    return collect_result_csv(summary_dir, "*.csv")


def collect_analysis_names() -> list[str]:
    """Discover all DEG/GSEA/TF analysis names that should become report blocks."""
    names: set[str] = set()
    for analysis_dir in dirs(TABLE_ROOT):
        if analysis_dir.name == "GSEA_summary":
            continue
        if (analysis_dir / "DEG").exists() or (analysis_dir / "GSEA").exists():
            names.add(analysis_dir.name)
    for type_dir_name in ("DEG", "deg"):
        type_dir = TF_SUMMARY_ROOT / type_dir_name
        if type_dir.exists():
            names.update(path.name for path in dirs(type_dir))
    return sorted(names, key=sort_analysis_name)


def collect_intersection_names() -> list[str]:
    """Discover all intersect schemes from 02号、GSEA and 09号TF整合 outputs."""
    names: set[str] = set()
    if INTERSECT_ROOT.exists():
        names.update(path.name for path in dirs(INTERSECT_ROOT))
    for type_dir_name in ("intersect", "INTERSECT"):
        type_dir = TF_SUMMARY_ROOT / type_dir_name
        if type_dir.exists():
            names.update(path.name for path in dirs(type_dir))
    for analysis_dir in dirs(TABLE_ROOT):
        if "_" in analysis_dir.name and (analysis_dir / "GSEA").exists():
            names.add(analysis_dir.name)
    return sorted(names, key=sort_analysis_name)


def append_deg_frames(lines: list[str], analysis: str) -> None:
    """Append 01号差异分析 result pages for one analysis design."""
    deg_dir = TABLE_ROOT / analysis / "DEG"
    if not deg_dir.exists():
        return

    lines += text_frame(
        f"{analysis}：差异表达分析阅读顺序",
        [
            f"本组以 {analysis} 为一个完整分析单元，连续展示 summary、显著差异基因表、传统火山图和 Top DEG 表达热图。",
            "summary 用于确认显著基因总量、上下调方向和阈值；显著基因表用于定位具体候选；火山图查看 logFC 与显著性分布；Top DEG 热图回到样本表达层面验证分组分离。",
            "这四页共同回答：该分析设计中 LRC 相对 BULK 是否形成稳定、可解释、可被后续富集承接的休眠样转录改变。",
        ],
    )

    summary_csv = find_result_csv(deg_dir, "summary")
    if summary_csv is not None:
        tex = make_preview_table(summary_csv)
        if tex is not None:
            lines += table_frame(
                f"{analysis} | DEG统计摘要",
                tex,
                [
                    "该 summary 是当前 LRC/BULK 比较的入口页。重点核对 Up、Down 与 Total_Significant_Genes，三者应满足 Total=Up+Down。",
                    "P_Value_Column、P_Value_Cutoff 与 LogFC_Cutoff 说明当前脚本采用的显著性判定标准；这些阈值也会被火山图和后续富集输入沿用。",
                    "Samples_Used 帮助判断样本量是否足够支撑该比较；如果显著基因数很少或方向严重失衡，需要结合后续图形谨慎解释。",
                ],
                source_csv=summary_csv,
            )

    significant_csv = find_result_csv(deg_dir, "significant_genes")
    if significant_csv is not None:
        tex = make_preview_table(significant_csv)
        if tex is not None:
            lines += table_frame(
                f"{analysis} | 显著差异基因表",
                tex,
                [
                    "该表展示通过阈值筛选后的显著 DEG 前 21 行。Symbol/Ensembl/Entrez 用于后续交集、GSEA leading-edge 回查和 TF 富集。",
                    "logFC 为 LRC 相对 BULK 的效应量；正值表示 LRC 方向上调，负值表示 BULK 方向更高。t、P.Value、adj.P.Val 用于判断统计强度。",
                    "阅读时优先关注效应量大、校正后 p 值稳定，并且在多个模型或交集方案中反复出现的候选基因。",
                ],
                source_csv=significant_csv,
            )

    volcano_fig = PLOT_ROOT / "volcano" / analysis / "png" / "volcano_plot.png"
    if volcano_fig.exists():
        lines += figure_frame(
            f"{analysis} | 传统火山图",
            volcano_fig,
            [
                "火山图横坐标为 logFC，纵坐标为显著性指标；红色 Sig_Up 与蓝色 Sig_Down 分别表示 LRC 方向上调和下调显著基因。",
                "该页用于快速判断当前分析的显著基因是否呈现清晰方向结构：红蓝点数量、左右分布和标注基因的位置都可提示该模型的差异强度。",
                "若被标注基因位于高显著性、高效应量区域，并且在交集、GSEA 或 TF 整合中反复出现，应优先进入后续机制验证。",
            ],
        )

    heatmap_fig = PLOT_ROOT / "gene_heatmap" / analysis / "png" / "gene_heatmap.png"
    if heatmap_fig.exists():
        lines += figure_frame(
            f"{analysis} | Top DEG表达热图",
            heatmap_fig,
            [
                "Top DEG 热图将统计结果重新投射到样本表达矩阵上，检查这些候选基因是否能把 LRC 与 BULK 样本分开。",
                "顶部条带标记 Group 与 Cell_Line，左侧方向条保留 Up/Down 信息。若 LRC 样本在热图中形成一致表达块，说明 DEG 不只是统计显著，也具有样本层面的状态区分能力。",
                "该图是连接 DEG 与后续机制分析的重要质量控制：分离越清晰，后续用这些基因做交集、TF 富集或候选验证越有解释力。",
            ],
            wide=True,
        )


def append_intersection_frames(lines: list[str], scheme: str) -> None:
    """Append 02号显著基因交集 result pages for one intersect scheme."""
    scheme_dir = INTERSECT_ROOT / scheme
    if not scheme_dir.exists():
        return

    lines += text_frame(
        f"{scheme}：交集方案阅读顺序",
        [
            f"本组以 {scheme} 为一个完整交集单元，先展示交集统计摘要，再展示交集基因注释列表，随后回查参与交集的各 DEG 结果，最后展示多组火山图。",
            "交集策略的目的不是取最大基因数，而是从多个 LRC/BULK 模型中筛出方向更稳定、可重复性更强的候选基因。",
            "阅读时要同时看交集规模、方向一致性和成员 DEG 统计量；只有在多个成员分析中方向和显著性都较稳定的基因，才更适合进入后续机制解释。",
        ],
    )

    summary_csv = find_result_csv(scheme_dir, "summary")
    if summary_csv is not None:
        tex = make_preview_table(summary_csv)
        if tex is not None:
            lines += table_frame(
                f"{scheme} | 交集统计摘要",
                tex,
                [
                    "该 summary 用于判断当前交集方案的严格程度。Total_Intersected_Genes 是交集基因数量；Common_Up/Common_Down 表示方向一致候选；Mixed_Direction 提示成员分析间方向不完全一致。",
                    "如果交集基因很少但方向一致，说明该组合更严格、更偏向稳定候选；如果 Mixed_Direction 较多，则需要回到成员 DEG 表逐个判断是否值得保留。",
                ],
                source_csv=summary_csv,
            )

    gene_list_csv = find_result_csv(scheme_dir, "gene_list")
    if gene_list_csv is not None:
        tex = make_preview_table(gene_list_csv)
        if tex is not None:
            lines += table_frame(
                f"{scheme} | 交集基因注释列表",
                tex,
                [
                    "该 gene_list 只保留交集基因的可追踪注释信息，适合作为后续 TF 富集、人工筛选、文献检索和实验验证的输入。",
                    "这里不混合展示 p 值和 logFC，是为了避免不同 DEG 设计的统计量被误读；统计证据会在后续成员 DEG 结果中分别回查。",
                ],
                source_csv=gene_list_csv,
            )

    member_deg_files = sorted(
        collect_result_csv(scheme_dir, "deg_results.csv"),
        key=lambda path: sort_analysis_name(normalize_result_csv_path(path).parent.name),
    )
    for csv_path in member_deg_files:
        flat_csv = normalize_result_csv_path(csv_path)
        rel_parts = flat_csv.relative_to(scheme_dir).parts
        member = rel_parts[0] if len(rel_parts) > 1 else scheme
        tex = make_preview_table(csv_path)
        if tex is None:
            continue
        lines += table_frame(
            f"{scheme} | {member}中的交集基因DEG结果",
            tex,
            [
                f"该表回查 {scheme} 交集基因在 {member} 原始 DEG 结果中的 logFC、t、P.Value 和 adj.P.Val。",
                "如果同一 Symbol 在不同成员分析中 logFC 方向一致，并且显著性稳定，说明它更可能代表跨模型休眠样转录程序，而不是单一细胞系背景噪音。",
            ],
            source_csv=csv_path,
        )

    multiple_volcano_fig = PLOT_ROOT / "multiple_volcano" / scheme / "png" / "multiple_volcano_plot.png"
    if multiple_volcano_fig.exists():
        lines += figure_frame(
            f"{scheme} | 多组火山图",
            multiple_volcano_fig,
            [
                "多组火山图把参与该交集方案的多个 LRC/BULK 分析并列展示，红色 Sig_Up 与蓝色 Sig_Down 对应各模型显著上/下调基因。",
                "该图用于判断不同模型中的显著 DEG 是否具有相似方向结构。若多个成员均显示相近的红蓝分布和强效候选，说明该交集方案更可能提炼到稳定生物信号。",
                "中心组名与每组坐标轴帮助快速定位参与比较的模型，后续 GSEA 和 TF 整合可继续沿用该交集方案进行机制解释。",
            ],
            wide=True,
        )


def append_gsea_frames(lines: list[str], analysis: str) -> None:
    """Append paired 06/07号GSEA dotplot and result-table pages for one scheme."""
    gsea_dir = TABLE_ROOT / analysis / "GSEA"
    if not gsea_dir.exists():
        return

    csv_files = sorted(collect_result_csv(gsea_dir, "gsea_result.csv"), key=sort_gsea_csv)
    if not csv_files:
        return

    lines += text_frame(
        f"{analysis}：GSEA通路富集阅读顺序",
        [
            "GSEA 使用全量 ranked genes，而不是只使用显著基因列表；rank metric 由 06 号脚本统一定义，当前用于判断 LRC/BULK 状态的整体通路偏移。",
            "每个 MSigDB 类别先展示 dotplot，再展示同一结果表。dotplot 用于快速看 top10 通路主题，结果表用于核对 NES、pvalue、p.adjust 和 qvalue。",
            "NES 为正通常表示 LRC 方向富集，NES 为负通常表示 BULK 方向富集。优先关注与炎症、免疫微环境、ECM remodeling、细胞周期、代谢重编程和 TF 靶集相关的通路。",
        ],
    )

    preview_map = make_preview_tables_parallel(csv_files)
    for csv_path in csv_files:
        flat_csv = normalize_result_csv_path(csv_path)
        geneset = flat_csv.relative_to(TABLE_ROOT / analysis / "GSEA").parts[0]
        fig = PLOT_ROOT / "GSEA" / analysis / geneset / "png" / "dotplot.png"
        if fig.exists():
            lines += figure_frame(
                f"{analysis} | GSEA dotplot | {gsea_display_label(geneset)}",
                fig,
                [
                    "该 dotplot 展示当前 analysis × MSigDB 类别中通过统一阈值后的 top10 通路，用于快速定位最强富集主题。",
                    "气泡大小和颜色同时表达富集规模和统计显著性；通路名称若与休眠、苏醒、免疫、炎症、ECM、细胞周期或 TF 靶集相关，应优先进入后续解释。",
                    "下一页为同一 GSEA 设计的统计表，可用 NES、p.adjust 与 qvalue 对图中通路进行定量核对。",
                ],
                wide=True,
            )
        compact_tex = preview_map.get(resolve_result_csv_path(csv_path))
        if compact_tex is None:
            continue
        lines += table_frame(
            f"{analysis} | GSEA结果表 | {gsea_display_label(geneset)}",
            compact_tex,
            [
                "该表是上一页 dotplot 对应的 GSEA 统计结果预览。ID 为通路条目；NES 判断富集方向；pvalue、p.adjust 和 qvalue 判断统计证据。",
                "本报告展示核心列以便汇报阅读；完整 CSV 中仍保留 clusterProfiler 官方输出的其他字段，可用于回查 leading-edge 基因和通路细节。",
                "当某通路同时具有明确 NES 方向和稳定校正显著性时，可作为休眠维持或休眠细胞苏醒机制的后续重点。",
            ],
            source_csv=csv_path,
        )


def append_tf_summary_frames(lines: list[str], input_type: str, analysis: str) -> None:
    """Append 09号TF整合 pages for one DEG/intersect input scheme."""
    csv_files = collect_tf_summary_csv(input_type, analysis)
    if not csv_files:
        return

    input_label = "DEG" if input_type.lower() == "deg" else "intersect"
    lines += text_frame(
        f"{analysis}：TF整合结果阅读顺序",
        [
            f"本组 TF 结果来自 {input_label}/{analysis} 输入方案。09 号脚本将六种 TF 方法的原始输出整理为 method_final 排名，并进一步计算多方法交集。",
            "先看 TF 交集方案摘要，判断不同方法组合的候选数量；再看六种单方法 final 排名，理解每种方法的独立证据；最后看 Top10 交集候选，筛选更稳健的上游调控因子。",
            "Source_Method_Count 越高、CheA3_Library_Count 越大，通常说明跨方法和外部 library 支持更充分；VIPER/CollecTRI 的方向信息可辅助判断 TF 活性倾向。",
        ],
    )

    summary_dir = existing_tf_summary_dir(input_type, analysis)
    if summary_dir is None:
        return
    csv_files = sorted(csv_files, key=sort_tf_csv)
    preview_map = make_preview_tables_parallel(csv_files)

    method_summary: Path | None = None
    method_results: list[Path] = []
    intersection_summary: Path | None = None
    scheme_candidates: dict[str, Path] = {}

    for csv_path in csv_files:
        flat_csv = normalize_result_csv_path(csv_path)
        parts = flat_csv.relative_to(summary_dir).parts
        if parts == ("method_final", "summary", "method_final_summary.csv"):
            method_summary = csv_path
        elif len(parts) >= 3 and parts[0] == "method_final" and parts[1] != "summary":
            method_results.append(csv_path)
        elif parts == ("intersections", "summary", "intersection_summary.csv"):
            intersection_summary = csv_path
        elif len(parts) >= 4 and parts[0] == "intersections" and parts[2] == "candidates":
            scheme_candidates[parts[1]] = csv_path

    if intersection_summary is not None:
        tex = preview_map.get(resolve_result_csv_path(intersection_summary))
        if tex is not None:
            lines += table_frame(
                f"{analysis} | TF交集方案摘要",
                tex,
                [
                    "该表汇总当前输入方案下全部 TF 交集组合。Intersection_Name 表示组合名；Number_Of_Methods 为参与方法数量；Intersected_TF_Count 为全量交集 TF 数。",
                    "该页用于判断哪些组合更严格、哪些组合保留候选更多。严格组合若仍有候选，通常更适合优先验证；宽松组合则适合提供补充线索。",
                ],
                source_csv=intersection_summary,
            )

    if method_summary is not None:
        tex = preview_map.get(resolve_result_csv_path(method_summary))
        if tex is not None:
            lines += table_frame(
                f"{analysis} | 六方法final结果摘要",
                tex,
                [
                    "该表展示六种 TF 富集/活性推断方法最终进入排序和交集分析的 TF 数量。",
                    "如果某方法 Final_TF_Count 很低，说明该方法在当前输入基因集或 ranked signature 上覆盖有限；解释交集结果时要考虑方法覆盖度。",
                ],
                source_csv=method_summary,
            )

    for csv_path in sorted(method_results, key=sort_tf_csv):
        flat_csv = normalize_result_csv_path(csv_path)
        method = flat_csv.parent.name
        tex = preview_map.get(resolve_result_csv_path(csv_path))
        if tex is None:
            continue
        lines += table_frame(
            f"{analysis} | 单方法TF final排序 | {tf_method_label(method)}",
            tex,
            [
                f"该表展示 {tf_method_label(method)} 方法整理后的 final TF 排名。不同方法的 Score、P_Value、NES 或 Direction 含义不同，应结合方法类型理解。",
                "Rank/TF 是最直接的候选定位字段；CheA3_Library_Count 与 CheA3_Integrated_TopRank 用于补充跨数据库可靠性判断。",
                "单方法排名用于保留方法特异线索，最终候选优先结合后续多方法交集表判断。",
            ],
            source_csv=csv_path,
        )

    ordered_schemes = sorted(
        set(scheme_candidates),
        key=lambda scheme: (TF_INTERSECTION_ORDER.get(scheme, 99), scheme),
    )
    candidate_entries: list[tuple[str, Path, Path]] = []
    for scheme in ordered_schemes:
        csv_path = scheme_candidates.get(scheme)
        if csv_path is None:
            continue
        tex = preview_map.get(resolve_result_csv_path(csv_path))
        if tex is not None:
            candidate_entries.append((scheme, tex, csv_path))

    candidate_note = [
        "Consensus_Rank 为综合排序；TF 为转录因子 symbol；Source_Method_Count 表示支持该 TF 的方法数。",
        "Source_Methods 显示来自哪些方法；CheA3_Library_Count 表示 ChEA3 多 library 支持数量，越高通常说明外部证据越丰富。",
        "优先关注多方法支持、ChEA3 library 较多，并且在 VIPER/CollecTRI 等活性方法中方向较清楚的 TF。",
    ]
    for i in range(0, len(candidate_entries), 2):
        first = candidate_entries[i]
        second = candidate_entries[i + 1] if i + 1 < len(candidate_entries) else None
        first_title = tf_intersection_label(first[0])
        first_items = [f"上表：{first_title}。"] + candidate_note
        if second is None:
            lines += table_frame(
                f"{analysis} | TF交集Top10候选 | {first_title}",
                first[1],
                first_items,
                source_csv=first[2],
            )
            continue
        second_title = tf_intersection_label(second[0])
        second_items = [f"下表：{second_title}。"] + candidate_note
        lines += stacked_two_table_frame(
            f"{analysis} | TF交集Top10候选 | {first_title} / {second_title}",
            first[1],
            first[2],
            first_items,
            second[1],
            second[2],
            second_items,
        )


def build_result_overview() -> str:
    """Build one compact page explaining the new analysis-centric report logic."""
    section_name = f"{DATASET_ID}/00_report_logic_overview.tex"
    lines = section_cover(
        "结果展示逻辑：按分析方案串联证据链",
        "从单个 DEG 方案或交集方案出发，连续查看 DEG、图形、GSEA 和 TF 证据",
        "本版报告不再严格按脚本编号线性展示，而是把每个差异分析方案、每个交集方案作为一个独立证据单元。这样更接近汇报时的阅读方式：先定义一个模型/交集，再连续判断它的基因、通路和上游调控因子。",
    )
    lines += text_frame(
        "本报告的主线结构",
        [
            "第一部分仍保留研究设计与 00 号样本结构质控，用于说明课题逻辑和数据基础。",
            "随后进入 DEG 分析方案主线：每个 analysis 依次展示 DEG summary、显著基因、传统火山图、Top DEG 热图、GSEA 图表和 TF 整合结果。",
            "再进入 intersect 交集方案主线：每个交集方案依次展示交集 summary、交集基因注释、成员 DEG 结果、多组火山图、GSEA 图表和 TF 整合结果。",
            "这种结构有助于在同一方案内完整追踪证据：从差异基因到通路，再到可能驱动休眠维持或苏醒的 TF。",
        ],
    )
    write_text(section_file(section_name), lines)
    return section_name


def build_analysis_scheme_sections() -> list[str]:
    """Build Beamer sections using DEG analysis names as the primary order."""
    section_names: list[str] = []
    for analysis in collect_analysis_names():
        section_name = f"{DATASET_ID}/10_analysis_scheme_{sanitize_name(analysis)}.tex"
        section_names.append(section_name)
        lines = section_cover(
            f"分析方案：{analysis}",
            "DEG → 火山图/热图 → GSEA → TF整合 的连续证据链",
            f"本节围绕 {analysis} 这一差异分析方案展开。先确认 LRC/BULK 差异基因是否稳定，再查看通路层面的方向性富集，最后整合多方法 TF 结果寻找可能驱动休眠样状态的上游调控因子。",
        )
        append_deg_frames(lines, analysis)
        append_gsea_frames(lines, analysis)
        append_tf_summary_frames(lines, "DEG", analysis)
        write_text(section_file(section_name), lines)
    return section_names


def build_intersection_scheme_sections() -> list[str]:
    """Build Beamer sections using intersect schemes as the primary order."""
    section_names: list[str] = []
    for scheme in collect_intersection_names():
        section_name = f"{DATASET_ID}/20_intersection_scheme_{sanitize_name(scheme)}.tex"
        section_names.append(section_name)
        lines = section_cover(
            f"交集方案：{scheme}",
            "交集基因 → 多组火山图 → GSEA → TF整合 的跨模型证据链",
            f"本节围绕 {scheme} 这一交集方案展开。它的核心目的，是从多个 LRC/BULK 模型中提炼方向更稳定、可重复性更高的候选基因，并进一步追踪这些候选对应的通路和 TF 证据。",
        )
        append_intersection_frames(lines, scheme)
        append_gsea_frames(lines, scheme)
        append_tf_summary_frames(lines, "intersect", scheme)
        write_text(section_file(section_name), lines)
    return section_names


def build_main(section_names: list[str]) -> None:
    theme_path = THEME_DIR.as_posix()
    lines = [
        "% Auto-built Beamer report for CRLM_dormancy / GSE114012",
        "% Source generator: scripts/functions/beamer_dataset_report_builder.py",
        r"\makeatletter",
        rf"\def\input@path{{{{{theme_path}/}}}}",
        r"\makeatother",
        "",
        r"\documentclass[8pt,aspectratio=169,xcolor=dvipsnames]{beamer}",
        r"\usetheme{SimplePlus}",
        "",
        r"\usepackage[UTF8]{ctex}",
        r"\usepackage{fontspec}",
        r"\usepackage{graphicx}",
        r"\usepackage{booktabs}",
        r"\usepackage{array}",
        r"\usepackage{adjustbox}",
        r"\usepackage{ragged2e}",
        r"\usepackage{xcolor}",
        r"\usepackage{newunicodechar}",
        r"\usepackage{bookmark}",
        r"\setsansfont{Helvetica}",
        r"\newunicodechar{→}{\ensuremath{\rightarrow}}",
        r"\newunicodechar{κ}{\ensuremath{\kappa}}",
        r"\setbeamertemplate{navigation symbols}{}",
        r"\setbeamersize{text margin left=1.5mm,text margin right=1.5mm}",
        r"\newcommand{\ReportBodyFont}{\fontsize{6.55pt}{8.15pt}\selectfont}",
        r"\setbeamerfont{normal text}{size*={6.55pt}{8.15pt}}",
        r"\setbeamerfont{frametitle}{size=\Large,series=\bfseries}",
        r"\AtBeginEnvironment{frame}{\ReportBodyFont}",
        r"\hfuzz=80pt",
        r"\vfuzz=80pt",
        r"\hbadness=10000",
        r"\vbadness=10000",
        "",
        rf"\newcommand{{\CRLMROOT}}{{{PROJECT_ROOT.as_posix()}}}",
        r"\newcommand{\PathText}[1]{\begingroup\tiny\ttfamily\detokenize{#1}\endgroup}",
        r"\newcommand{\MissingBox}[1]{\fbox{\begin{minipage}[c][0.62\textheight][c]{0.92\linewidth}\centering\ttfamily\scriptsize Missing output\\[0.4em]\detokenize{#1}\end{minipage}}}",
        r"\newsavebox{\ResultImageBox}",
        r"\newlength{\ResultImageTextWidth}",
        "",
        r"\newcommand{\SectionCover}[3]{%",
        r"  \section{#1}%",
        r"  \begin{frame}[plain,t]%",
        r"    \vspace{1.2mm}%",
        r"    {\Large\bfseries #1}\par\vspace{0.45em}%",
        r"    {\normalsize #2}\par\vspace{0.8em}%",
        r"    \noindent\begin{minipage}{0.96\textwidth}\RaggedRight\ReportBodyFont #3\end{minipage}%",
        r"  \end{frame}%",
        r"}",
        "",
        r"\newcommand{\TextFrame}[2]{%",
        r"  \begin{frame}[t]{#1}%",
        r"    \RaggedRight\ReportBodyFont #2%",
        r"  \end{frame}%",
        r"}",
        "",
        r"\newcommand{\ResultFigureFrame}[3]{%",
        r"  \begin{frame}[plain,t]%",
        r"    {\normalsize\bfseries #1}\par\vspace{0.06em}%",
        r"    \IfFileExists{\CRLMROOT/\detokenize{#2}}{%",
        r"      \sbox{\ResultImageBox}{\includegraphics[width=0.69\textwidth,height=0.915\paperheight,keepaspectratio]{\CRLMROOT/\detokenize{#2}}}%",
        r"      \setlength{\ResultImageTextWidth}{\dimexpr\textwidth-\wd\ResultImageBox-0.6mm\relax}%",
        r"      \noindent\parbox[t]{\wd\ResultImageBox}{\vspace{0pt}\noindent\usebox{\ResultImageBox}}%",
        r"      \parbox[t]{\ResultImageTextWidth}{\vspace{0pt}\RaggedRight\ReportBodyFont\setlength{\rightskip}{0pt plus 1fil}#3}%",
        r"    }{\MissingBox{#2}}%",
        r"  \end{frame}%",
        r"}",
        "",
        r"\newcommand{\ResultWideFigureFrame}[3]{%",
        r"  \begin{frame}[plain,t]%",
        r"    {\normalsize\bfseries #1}\par\vspace{0.06em}%",
        r"    \IfFileExists{\CRLMROOT/\detokenize{#2}}{%",
        r"      \sbox{\ResultImageBox}{\includegraphics[width=0.735\textwidth,height=0.915\paperheight,keepaspectratio]{\CRLMROOT/\detokenize{#2}}}%",
        r"      \setlength{\ResultImageTextWidth}{\dimexpr\textwidth-\wd\ResultImageBox-0.6mm\relax}%",
        r"      \noindent\parbox[t]{\wd\ResultImageBox}{\vspace{0pt}\noindent\usebox{\ResultImageBox}}%",
        r"      \parbox[t]{\ResultImageTextWidth}{\vspace{0pt}\RaggedRight\ReportBodyFont\setlength{\rightskip}{0pt plus 1fil}#3}%",
        r"    }{\MissingBox{#2}}%",
        r"  \end{frame}%",
        r"}",
        "",
        r"\newcommand{\ResultTableFrame}[4]{%",
        r"  \begin{frame}[plain,t]{#1}%",
        r"    \vspace{-0.28em}%",
        r"    \begin{minipage}[t]{\textwidth}\RaggedRight\ReportBodyFont #4\end{minipage}%",
        r"    \par\vspace{0.22em}%",
        r"    \IfFileExists{\CRLMROOT/\detokenize{#3}}{%",
        r"      \begingroup\centering%",
        r"      \begin{adjustbox}{max width=0.948\textwidth,max totalheight=0.755\textheight,keepaspectratio}%",
        r"        \input{\CRLMROOT/\detokenize{#3}}%",
        r"      \end{adjustbox}%",
        r"      \par\endgroup%",
        r"    }{\MissingBox{#3}}%",
        r"  \end{frame}%",
        r"}",
        "",
        r"\newcommand{\ResultTwoTableFrame}[7]{%",
        r"  \begin{frame}[plain,t]{#1}%",
        r"    \vspace{-0.28em}%",
        r"    \noindent\begin{minipage}[t]{0.498\textwidth}\vspace{0pt}%",
        r"      \RaggedRight\ReportBodyFont #4%",
        r"      \par\vspace{0.18em}%",
        r"      \IfFileExists{\CRLMROOT/\detokenize{#3}}{%",
        r"        \begingroup\centering%",
        r"        \begin{adjustbox}{max width=0.955\linewidth,max totalheight=0.632\textheight,keepaspectratio}%",
        r"          \input{\CRLMROOT/\detokenize{#3}}%",
        r"        \end{adjustbox}%",
        r"        \par\endgroup%",
        r"      }{\MissingBox{#3}}%",
        r"    \end{minipage}%",
        r"    \hfill%",
        r"    \begin{minipage}[t]{0.498\textwidth}\vspace{0pt}%",
        r"      \RaggedRight\ReportBodyFont #7%",
        r"      \par\vspace{0.18em}%",
        r"      \IfFileExists{\CRLMROOT/\detokenize{#6}}{%",
        r"        \begingroup\centering%",
        r"        \begin{adjustbox}{max width=0.955\linewidth,max totalheight=0.632\textheight,keepaspectratio}%",
        r"          \input{\CRLMROOT/\detokenize{#6}}%",
        r"        \end{adjustbox}%",
        r"        \par\endgroup%",
        r"      }{\MissingBox{#6}}%",
        r"    \end{minipage}%",
        r"  \end{frame}%",
        r"}",
        "",
        r"\newcommand{\ResultStackedTwoTableFrame}[7]{%",
        r"  \begin{frame}[plain,t]{#1}%",
        r"    \vspace{-0.28em}%",
        r"    \begin{minipage}[t]{\textwidth}\vspace{0pt}%",
        r"      \RaggedRight\ReportBodyFont #4%",
        r"      \par\vspace{0.16em}%",
        r"      \IfFileExists{\CRLMROOT/\detokenize{#3}}{%",
        r"        \begingroup\centering%",
        r"        \begin{adjustbox}{max width=0.948\textwidth,max totalheight=0.302\textheight,keepaspectratio}%",
        r"          \input{\CRLMROOT/\detokenize{#3}}%",
        r"        \end{adjustbox}%",
        r"        \par\endgroup%",
        r"      }{\MissingBox{#3}}%",
        r"    \end{minipage}%",
        r"    \par\vspace{0.18em}%",
        r"    \begin{minipage}[t]{\textwidth}\vspace{0pt}%",
        r"      \RaggedRight\ReportBodyFont #7%",
        r"      \par\vspace{0.16em}%",
        r"      \IfFileExists{\CRLMROOT/\detokenize{#6}}{%",
        r"        \begingroup\centering%",
        r"        \begin{adjustbox}{max width=0.948\textwidth,max totalheight=0.302\textheight,keepaspectratio}%",
        r"          \input{\CRLMROOT/\detokenize{#6}}%",
        r"        \end{adjustbox}%",
        r"        \par\endgroup%",
        r"      }{\MissingBox{#6}}%",
        r"    \end{minipage}%",
        r"  \end{frame}%",
        r"}",
        "",
        r"\title[CRC Dormancy]{结直肠癌休眠细胞复苏的转录组证据链}",
        r"\subtitle{GSE114012 LRC/BULK 模型的可复现分析进展}",
        r"\author{CRLM dormancy project}",
        r"\institute{Internal exploratory report}",
        r"\date{\today}",
        "",
        r"\begin{document}",
        r"\begin{frame}\titlepage\end{frame}",
    ]
    lines.extend([
        rf"\input{{\detokenize{{{(SECTION_ROOT / name).as_posix()}}}}}"
        for name in section_names
    ])
    lines.extend([
        r"\end{document}",
        "",
    ])
    write_text(MAIN_TEX, lines)


def clean_generated_sections() -> None:
    # 保留已有文件，依靠write_text的内容比较实现“未变化不重写”。
    # 这样后续latexmk可以复用temporary/beamer中的aux状态，避免每次从零编译。
    SECTION_ROOT.mkdir(parents=True, exist_ok=True)
    PROJECT_SECTION_ROOT.mkdir(parents=True, exist_ok=True)
    DATASET_SECTION_ROOT.mkdir(parents=True, exist_ok=True)
    GENERATED_TABLE_ROOT.mkdir(parents=True, exist_ok=True)


def remove_stale_generated_files(section_names: list[str]) -> None:
    """Remove generated Beamer fragments that are no longer referenced.

    This keeps the report dynamic: if a future GSEA run only includes a subset
    of MSigDB collections, old table/dotplot section files from previous runs
    will not remain in the source tree.
    """
    expected_sections = {(SECTION_ROOT / name).resolve() for name in section_names}
    for root in (PROJECT_SECTION_ROOT, DATASET_SECTION_ROOT):
        for tex_path in root.glob("*.tex"):
            if tex_path.resolve() not in expected_sections:
                tex_path.unlink()

    for tex_path in GENERATED_TABLE_ROOT.glob("*.tex"):
        if tex_path.resolve() not in GENERATED_PREVIEW_TABLES:
            tex_path.unlink()


def generate_sources() -> list[str]:
    BEAMER_ROOT.mkdir(parents=True, exist_ok=True)
    clean_generated_sections()
    section_names: list[str] = []
    section_names.append(build_project_design())
    section_names.append(build_00())
    section_names.append(build_result_overview())
    section_names.extend(build_analysis_scheme_sections())
    section_names.extend(build_intersection_scheme_sections())
    build_main(section_names)
    remove_stale_generated_files(section_names)
    return section_names


def compile_report() -> None:
    BUILD_ROOT.mkdir(parents=True, exist_ok=True)
    REPORT_ROOT.mkdir(parents=True, exist_ok=True)
    command = [
        "latexmk",
        "-xelatex",
        "-silent",
        "-synctex=1",
        "-interaction=nonstopmode",
        "-file-line-error",
        f"-outdir={BUILD_ROOT.as_posix()}",
        f"-auxdir={BUILD_ROOT.as_posix()}",
        "-shell-escape",
        MAIN_TEX.as_posix(),
    ]
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
    built_pdf = BUILD_ROOT / "beamer_report.pdf"
    if not built_pdf.exists():
        raise FileNotFoundError(f"Compiled PDF was not found: {built_pdf}")
    shutil.copy2(built_pdf, FINAL_PDF)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build GSE114012 Beamer sources.")
    parser.add_argument("--compile", action="store_true", help="Compile the generated Beamer report.")
    args = parser.parse_args()

    sections = generate_sources()
    print(f"Generated main tex: {rel(MAIN_TEX)}")
    print(f"Generated section files: {len(sections)}")
    if args.compile:
        compile_report()
        print(f"Compiled PDF: {rel(FINAL_PDF)}")


if __name__ == "__main__":
    main()
