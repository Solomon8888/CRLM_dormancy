#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build modular Beamer report sources for dataset-level CRLM analyses.

This script currently scans the analysis outputs under results/ngs/GSE114012
and regenerates a modular Beamer source tree:

  scripts/beamer/beamer_report.tex
  scripts/beamer/sections/*.tex

When called with --compile, it also compiles the deck with latexmk using:

  temporary/beamer            for intermediate build files
  results/reports/beamer      for the final PDF

The text blocks are intentionally stable and dataset-specific. Figures and
tables are discovered from the current result paths, so rerunning the analysis
scripts followed by this builder refreshes the report. The constants near the
top can be extended later for additional datasets with the same output layout.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
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

TABLE_PREVIEW_ROWS = 21
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
    "all_genes": 3,
    "deg_results": 4,
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
        "GeneID",
        "Symbol",
        "Ensembl",
        "Entrez",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B",
    ],
    "intersect_summary": [
        "Selected_Analyses",
        "DEG_Result_Analyses",
        "Total_Intersected_Genes",
        "Common_Up",
        "Common_Down",
        "Mixed_Direction",
    ],
    "intersect_gene_list": [
        "GeneID",
        "Symbol",
        "Ensembl",
        "Entrez",
    ],
    "intersect_deg_results": [
        "GeneID",
        "Symbol",
        "Ensembl",
        "Entrez",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B",
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
        "Description",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust",
        "qvalue",
        "rank",
        "leading_edge",
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
        "Required_Methods",
        "Source_Method_Count",
        "Source_Methods",
        "Mean_Selected_Rank",
        "Best_Selected_Rank",
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
    if "DEG" in parts and stem in {"all_genes", "significant_genes"}:
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


def format_table_cell(value: object, max_chars: int = 90) -> str:
    if value is None:
        return "--"
    value = str(value)
    if value == "" or value.lower() == "nan":
        return "--"
    try:
        number = float(value)
        if abs(number) >= 1000 or (abs(number) > 0 and abs(number) < 0.001):
            return f"{number:.3e}"
        return f"{number:.4g}"
    except ValueError:
        return shorten_text(value, max_chars=max_chars)


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
    col_spec = "l" * len(headers)
    lines = [
        r"\begingroup",
        r"\tiny",
        r"\setlength{\tabcolsep}{2pt}",
        r"\renewcommand{\arraystretch}{1.28}",
        rf"\begin{{tabular}}{{@{{}}{col_spec}@{{}}}}",
        r"\toprule",
        " & ".join(tex_escape(header) for header in headers) + r" \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(" & ".join(tex_escape(value) for value in row) + r" \\")
    lines.extend([
        r"\bottomrule",
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
        r"    \setlength{\itemsep}{0.28em}",
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
        r"      \setlength{\itemsep}{0.20em}",
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
        r"      \setlength{\itemsep}{0.16em}",
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
        r"      \setlength{\itemsep}{0.16em}",
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
        r"      \setlength{\itemsep}{0.10em}",
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
        r"      \setlength{\itemsep}{0.10em}",
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
        r"      \setlength{\itemsep}{0.20em}",
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
        "all_genes": "全量差异分析排序表",
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
    rows = [
        ("00", "00_sample_clustering_heatmap.R", "确认 LRC/BULK 样本结构是否支持后续比较"),
        ("01", "01_limma_differential_expression.R", "定义休眠样 LRC 相对 BULK 的差异转录程序"),
        ("02", "02_intersect_significant_genes.R", "提炼跨细胞系反复出现的稳定候选基因"),
        ("04", "04_multiple_volcano_plot.R", "并列比较多模型 DEG 方向与效应量分布"),
        ("05", "05_top_deg_gene_heatmap.R", "检验 Top DEG 是否在样本层面分离 LRC 与 BULK"),
        ("06", "06_gsea_analysis.R", "基于全量排序基因运行 MSigDB GSEA"),
        ("07", "07_gsea_plot.R", "按同一 GSEA 设计配对展示 dotplot 与结果表"),
        ("08", "08_tf_enrichment_analysis.R", "用六类 TF 方法生成原始调控因子证据"),
        ("09", "09_integrate_tf_enrichment_results.R", "整合 TF 方法 final 结果并计算候选交集"),
    ]
    lines: list[str] = []
    lines += section_cover(
        "研究设计：从术后复发到休眠细胞复苏",
        "以 GSE114012 的 LRC/BULK 模型建立结直肠癌休眠样转录组证据链",
        "本报告不是简单罗列脚本输出，而是把每一步结果放回“休眠癌细胞如何被重新激活”的课题假说中：先确认模型，再定义差异程序，再提炼稳定基因、通路和 TF 调控候选。",
    )
    lines += text_frame(
        "临床问题：远期复发可能来自休眠癌细胞的重新苏醒",
        [
            "结直肠癌患者在根治性治疗后仍可能多年甚至十余年后出现远处复发，这提示部分肿瘤细胞可能早期已经播散，并在远处器官或微环境中长期维持低增殖状态。",
            "本课题的核心问题不是“癌细胞是否已经播散”，而是“什么信号使休眠或低增殖癌细胞重新进入细胞周期并形成可检测转移灶”。",
            "因此，维持播散癌细胞休眠或阻断复苏信号，可能成为预防复发和转移的治疗策略；本报告当前阶段聚焦转录组层面的候选基因、通路和 TF 证据。",
        ],
    )
    lines += text_frame(
        "GSE114012 在本课题中的定位",
        [
            "GSE114012 的实验主题是 colorectal cancer spheroids 中 dormant-like cells 的识别；当前项目将其中 CFSE label-retaining cells 视作 LRC 休眠样/低增殖群体，将 cycling cells 视作 BULK 对照群体。",
            "该数据覆盖 DLD1、HCT15、HT55、SW948、RKO、SW48 六种 CRC 细胞系，并保留多重复样本，因此既可以做整体 LRC vs BULK，也可以做单细胞系或组合模型比较。",
            "本报告基于整理后的 SummarizedExperiment 对象、临床样本表和统一结果目录，依次展示样本结构、差异表达、交集、表达热图、GSEA 和 TF 整合。",
            "所有结果均按固定路径动态扫描生成 Beamer；后续重跑 R 脚本后，只需重建报告即可展示最新版本。",
        ],
    )
    lines += text_frame(
        "当前分析路线：从模型确认到候选调控轴",
        [
            "第一层证据是样本结构：TPM 相关性热图用于判断 LRC/BULK 内部是否存在稳定表达模式，以及是否存在明显离群样本。",
            "第二层证据是差异转录程序：DEG summary、显著基因和 Top DEG 热图用于定义 LRC 相对 BULK 的核心表达改变。",
            "第三层证据是稳定性：交集分析用于从多个细胞系或组合模型中提炼反复出现的候选基因，降低单模型偏差。",
            "第四层证据是机制解释：GSEA 提供通路层面的方向性，TF 富集/活性推断提供潜在上游调控因子，二者共同服务于后续调控轴和药物重定位。",
        ],
    )
    lines += text_frame(
        "核心科学假说",
        [
            "结直肠癌复发并不一定只是残留癌细胞持续增殖的结果，也可能来自长期休眠癌细胞在特定微环境信号刺激下重新苏醒。",
            "因此，阻断复苏信号或维持休眠状态，可能比单纯追求杀伤全部残留癌细胞更接近临床预防复发的长期策略。",
            "本项目把“休眠维持”和“苏醒启动”拆解为三个层面：癌细胞内在转录程序、外部通路/微环境刺激、以及可被药物干预的关键调控节点。",
        ],
    )
    lines += text_frame(
        "总体研究设计：模型、评分与状态识别",
        [
            "步骤一：系统收集 colorectal cancer dormancy、recurrence、minimal residual disease、liver metastasis relapse 相关数据集，优先纳入原发癌、转移癌、复发癌、治疗后残留癌细胞和休眠样癌细胞模型。",
            "步骤二：基于文献 dormancy markers 构建 dormancy score，并结合 proliferation score 区分低增殖休眠样状态与增殖活跃状态。",
            "步骤三：在单细胞或分选群体数据中识别休眠型、过渡型、苏醒型和增殖型细胞，核心比较设定为休眠型 vs 苏醒型。",
        ],
    )
    lines += text_frame(
        "总体研究设计：TF、轨迹与通路",
        [
            "步骤四：通过 SCENIC、DoRothEA、pySCENIC、ChIP-X enrichment 等方法寻找驱动苏醒的转录因子，重点关注多方法重复支持的 TF。",
            "步骤五：利用 Monocle3、Slingshot、CytoTRACE 构建休眠到苏醒再到增殖的 pseudotime 轨迹，沿轨迹识别动态 TF 和动态通路。",
            "步骤六：应用 GSVA/GSEA 解析苏醒相关通路，重点观察炎症、ECM remodeling、TGFβ/IL6、NFκB、WNT/β-catenin、YAP/STAT3 等通路是否被激活。",
        ],
    )
    lines += text_frame(
        "总体研究设计：微环境、调控轴与药物重定位",
        [
            "步骤七：分析免疫微环境重塑，重点关注 M2 macrophage、TREM2/APOE macrophage、CAF、Treg 和 exhausted CD8 是否伴随苏醒型癌细胞增加。",
            "步骤八：用 CellChat/NicheNet 解析癌细胞与免疫/基质细胞互作，提炼 CAF-derived TGFβ/IL6 → STAT3/YAP、Macrophage TNF/IL1β → NFκB、WNT niche → β-catenin/TCF 等候选轴。",
            "步骤九：将苏醒型上调基因、关键 TF 靶基因和 leading-edge genes 输入 CMap/LINCS/DGIdb/Enrichr Drug Signatures，筛选可逆转苏醒 signature、维持休眠状态的候选药物。",
            "步骤十：形成可验证假说：阻断关键 TF/通路可使残留癌细胞长期维持休眠，从而降低结直肠癌复发和远处转移风险。",
        ],
    )
    lines += compact_table_frame("脚本编号与章节组织", rows)
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
        "01. LRC/BULK 差异表达：定义休眠样转录程序",
        "从每个分析设计中提取 LRC 相对 BULK 的 DEG、方向和全量排序",
        "01 号结果是后续所有基因集、通路和 TF 推断的核心输入。这里的目标不是只得到一张 DEG 表，而是为不同细胞系背景下的休眠样状态建立可比较的 ranked gene signature。",
    )
    lines += text_frame(
        "差异表达结果如何进入后续证据链",
        [
            "summary 用于快速判断每个模型中 LRC/BULK 差异强度：上调、下调和总显著基因数决定后续候选空间大小。",
            "significant_genes 是 02 号交集、05 号 Top DEG 热图和 08 号 TF 富集的主要输入，用于寻找达阈值的候选基因。",
            "all_genes 保留全量排序信息，是 06/07 号 GSEA 的输入；GSEA 不能只依赖显著基因列表，而需要完整 ranked gene list。",
            "ALL 代表整体休眠样趋势，单细胞系与组合分析则帮助判断候选信号是否跨遗传背景可重复。",
        ],
    )
    table_records: list[tuple[str, str, Path, str]] = []
    for analysis_dir in dirs(TABLE_ROOT):
        deg_dir = analysis_dir / "DEG"
        if not deg_dir.exists():
            continue
        analysis = analysis_dir.name
        for stem, note in [
            ("summary", "该表是本分析设计的 DEG 统计总览。Up/Down/Total_Significant_Genes 用于判断 LRC 相对 BULK 的转录变化规模；阈值列用于确认所有结果采用一致判定标准。"),
            ("significant_genes", "该表展示通过阈值的显著 DEG。logFC 表示 LRC 相对 BULK 的方向和效应量；P.Value/adj.P.Val/B 等统计量用于评估差异证据强弱。"),
            ("all_genes", "该表展示全量基因排序结果。即使某些基因未达到显著阈值，它们仍参与 GSEA ranked-list 计算，因此该表决定通路富集方向和强度。"),
        ]:
            csv_path = find_result_csv(deg_dir, stem)
            if csv_path is not None:
                table_records.append((analysis, stem, csv_path, note))
    table_records.sort(key=lambda record: (RESULT_STEM_ORDER.get(record[1], 50), sort_analysis_name(record[0])))
    for analysis, stem, csv_path, note in table_records:
        tex = make_preview_table(csv_path)
        if tex is None:
            continue
        lines += table_frame(
            deg_table_title(analysis, stem),
            tex,
            [note, "表格为前 21 行汇报预览；完整 CSV 保存在当前结果目录。"],
            source_csv=csv_path,
        )
    write_text(section_file(section_name), lines)
    return section_name


def build_02() -> str:
    section_name = f"{DATASET_ID}/02_intersect_significant_genes.tex"
    lines = section_cover(
        "02. 跨模型稳定基因：从 DEG 中提炼可重复候选",
        "对多个 LRC/BULK 分析设计的显著基因取交集，并回看每个成员分析的 DEG 证据",
        "单个细胞系的 DEG 可能包含模型特异噪音；交集分析用于寻找在多个 CRC 背景下反复出现的休眠样候选基因，为后续 TF 富集、通路解释和实验验证提供更稳的基因集合。",
    )
    lines += text_frame(
        "如何阅读交集结果",
        [
            "每个交集方案对应一组分析设计，例如多个细胞系组合或两个模型之间的稳定候选提取。",
            "summary 先回答交集规模：有多少共同基因、共同上调、共同下调，以及是否存在方向不一致的 mixed-direction 基因。",
            "gene_list 只保留 Symbol/Ensembl/Entrez 等注释信息，适合作为后续富集或候选清单输入。",
            "成员分析下的 deg_results 用于回看这些交集基因在每个原始 DEG 结果中的 logFC、p 值和方向，防止只看交集名单而忽略统计证据。",
        ],
    )
    for scheme_dir in sorted(dirs(INTERSECT_ROOT), key=lambda path: sort_analysis_name(path.name)):
        scheme = scheme_dir.name
        for csv_path in sorted(collect_result_csv(scheme_dir, "*.csv"), key=sort_result_csv_by_stem):
            flat_csv = normalize_result_csv_path(csv_path)
            rel_parts = flat_csv.relative_to(scheme_dir).parts
            if rel_parts[-1] == "summary.csv":
                note = "该 summary 用于判断当前交集方案是否足够严格且仍保留可解释候选。Total_Intersected_Genes 是交集基因总数；Common_Up/Common_Down 表示方向一致的候选；Mixed_Direction 提示需要谨慎解释。"
            elif rel_parts[-1] == "gene_list.csv":
                note = "该 gene_list 是交集后的候选基因注释表，不包含 p 值和 logFC，适合作为后续通路富集、TF 富集或人工筛选的输入清单。"
            else:
                note = "该表把交集基因映射回对应 DEG 结果。logFC 用于看 LRC 相对 BULK 的方向；P.Value/adj.P.Val 用于看统计显著性；多个成员结果方向一致时，候选可信度更高。"
            tex = make_preview_table(csv_path)
            if tex is None:
                continue
            lines += table_frame(intersect_table_title(scheme, rel_parts), tex, [note], source_csv=csv_path)
    write_text(section_file(section_name), lines)
    return section_name


def build_03() -> str:
    section_name = f"{DATASET_ID}/03_volcano_plot.tex"
    lines = section_cover(
        "03. 传统火山图",
        "单个 DEG 设计的显著性与效应量",
        "03 号脚本批量读取 DEG all_genes，按 01 号脚本一致阈值绘制传统火山图。",
    )
    lines += text_frame(
        "03 号脚本逻辑",
        [
            "横坐标 logFC 左右对称，保证 0 点居中；纵坐标为显著性度量，红色为 Sig_Up，蓝色为 Sig_Down，灰色为 Not_Sig。",
            "若脚本头部配置指定基因，则标注指定基因；否则自动标注 top 上调和 top 下调基因。",
            "火山图用于快速判断每套分析中显著基因的数量、方向和效应量分布。",
        ],
    )
    figures = sorted(
        (PLOT_ROOT / "volcano").rglob("volcano_plot.png"),
        key=lambda path: sort_analysis_name(first_plot_level("volcano", path)),
    )
    for fig in figures:
        analysis = first_plot_level("volcano", fig)
        lines += figure_frame(
            image_result_title("03", analysis, "传统火山图"),
            fig,
            [
                "红蓝点分别表示通过阈值的上调与下调基因，Not_Sig 用灰色展示。",
                "该图帮助判断 LRC 相对 BULK 的差异幅度是否集中在少数强效基因，或呈现广泛转录重塑。",
            ],
        )
    write_text(section_file(section_name), lines)
    return section_name


def build_04() -> str:
    section_name = f"{DATASET_ID}/04_multiple_volcano_plot.tex"
    lines = section_cover(
        "04. 多模型 DEG 方向总览：哪些模型共享强差异信号",
        "用多组火山图并列观察不同 LRC/BULK 分析设计的显著上调与下调基因",
        "传统单图火山图不再纳入本报告；这里保留多组火山图，是因为它更适合比较多个细胞系或组合模型中的 DEG 方向一致性和效应量分布。",
    )
    lines += text_frame(
        "多组火山图在本课题中的用途",
        [
            "每个小面板对应一个 DEG 分析设计；红色表示 LRC 上调，蓝色表示 LRC 下调，灰色 NS 点不展示，从而突出真正进入候选空间的显著基因。",
            "纵坐标保留真实 logFC，因此可比较不同模型中显著基因的效应量大小；横向布局用于快速观察不同模型的显著基因密度。",
            "该图回答的是“哪些细胞系或组合模型中，休眠样 LRC 的转录改变更强、更集中、更方向一致”。",
            "后续如果某些 TF 或通路只在某个模型强烈出现，应回到本图判断它是否来自模型特异差异，而不是普遍休眠样程序。",
        ],
    )
    figures = sorted(
        (PLOT_ROOT / "multiple_volcano").rglob("multiple_volcano_plot.png"),
        key=lambda path: sort_analysis_name(first_plot_level("multiple_volcano", path)),
    )
    for fig in figures:
        analysis = first_plot_level("multiple_volcano", fig)
        lines += figure_frame(
            image_result_title("04", analysis, "多组火山图"),
            fig,
            [
                "红/蓝点分别代表在该分析设计中通过 DEG 阈值的 LRC 上调/下调基因；点的位置保留真实 logFC 信息。",
                "组名色块用于区分不同模型，目的是把多个 DEG 设计放在同一视觉坐标中比较，而不是单独解读某一个模型。",
                "若多个模型中显著基因数量、方向或效应量分布相似，说明该转录程序更可能代表共同休眠样状态；若差异很大，则提示细胞系背景可能影响 LRC 程序。",
            ],
            wide=True,
        )
    write_text(section_file(section_name), lines)
    return section_name


def build_05() -> str:
    section_name = f"{DATASET_ID}/05_top_deg_gene_heatmap.tex"
    lines = section_cover(
        "05. Top DEG 表达热图：候选基因能否分离 LRC 与 BULK",
        "用 Top50 DEG 的样本表达模式验证差异基因集的状态区分能力",
        "DEG 表只是统计结果；热图用于回到样本层面检查 Top DEG 是否真正形成 LRC/BULK 的表达分离。如果无法分离，后续通路和 TF 解读都需要降低确信度。",
    )
    lines += text_frame(
        "Top DEG 热图如何帮助解释 DEG",
        [
            "每套分析选择 Top50 显著差异基因，按样本的 log2(TPM+1) 表达计算 row z-score，突出同一基因在不同样本间的相对高低。",
            "实验组统一标记为 LRC，对照组标记为 BULK，并把 LRC 放在左侧、BULK 放在右侧，使读图方向与研究问题一致。",
            "顶部注释条显示 Group 和 Cell_Line；左侧方向条保留该基因在 DEG 中的 Up/Down 信息。",
            "如果 Top DEG 能清楚区分 LRC 与 BULK，说明这些候选基因不仅统计显著，也能在样本层面稳定表达分离。",
        ],
    )
    figures = sorted(
        (PLOT_ROOT / "gene_heatmap").rglob("gene_heatmap.png"),
        key=lambda path: sort_analysis_name(first_plot_level("gene_heatmap", path)),
    )
    for fig in figures:
        analysis = first_plot_level("gene_heatmap", fig)
        lines += figure_frame(
            image_result_title("05", analysis, "Top DEG表达热图"),
            fig,
            [
                "每行是一个 Top DEG，每列是一个样本；颜色代表该基因在所有样本中的标准化表达高低，而不是绝对 TPM。",
                "如果 LRC 样本在多个 Top DEG 上呈现一致高/低表达，并与 BULK 明显分离，说明该分析设计下的休眠样转录程序较稳定。",
                "右侧基因名和顶部样本注释可用于挑选后续验证对象，特别是同时出现在 02 号交集、06/07 号 leading-edge 或 09 号 TF 证据链中的基因。",
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
            "当前报告只保留 Hallmark、BioCarta、KEGG、Reactome、WikiPathways、TFT、GO、HPO、C6、ImmuneSigDB 等与机制解释相关的 MSigDB 类别。",
            "06 号表格给出每个 analysis × geneset 的运行规模、正/负 NES 通路数量和单通路图数量；07 号随后对同一设计进行 dotplot 与结果表配对解读。",
            "正 NES 通常提示 LRC 方向富集，负 NES 通常提示 BULK 方向富集；这为“休眠维持”或“苏醒启动”相关通路筛选提供方向性。",
        ],
    )
    summary_csv = find_result_csv(RESULT_ROOT / "tables" / "GSEA_summary", "summary")
    alt_summary_csv = find_result_csv(RESULT_ROOT / "GSEA" / "summary", "summary")
    if summary_csv is not None:
        summary_tex = make_preview_table(summary_csv)
        if summary_tex is not None:
            lines += table_frame("06. GSEA运行摘要 | 全部analysis与基因集类别", summary_tex, ["该 summary 展示每个 analysis × geneset 的 GSEA 运行概况。Ranked_Genes 为进入排序的基因数；GSEA_Terms 为通过当前显著性阈值的通路数；Positive_NES/Negative_NES 表示 LRC/BULK 方向富集数量。"], source_csv=summary_csv)
    elif alt_summary_csv is not None:
        summary_tex = make_preview_table(alt_summary_csv)
        if summary_tex is not None:
            lines += table_frame("06. GSEA运行摘要 | 全部analysis与基因集类别", summary_tex, ["该 summary 展示每个 analysis × geneset 的 GSEA 运行概况。Ranked_Genes 为进入排序的基因数；GSEA_Terms 为通过当前显著性阈值的通路数；Positive_NES/Negative_NES 表示 LRC/BULK 方向富集数量。"], source_csv=alt_summary_csv)
    write_text(section_file(section_names[0]), lines)

    return section_names


def build_07() -> list[str]:
    section_names = [f"{DATASET_ID}/07_gsea_plotting.tex"]
    lines = section_cover(
        "07. GSEA 图表配对解读：先看通路图，再看统计表",
        "每个 analysis/geneset 的 dotplot 与同一结果表连续展示",
        "本章节把每个 GSEA 设计的图片和表格放在一起阅读。dotplot 用于快速识别 top 通路，紧随其后的表格用于核对 NES、p.adjust、qvalue、rank 和 leading_edge 等统计字段。",
    )
    lines += text_frame(
        "GSEA 图表配对的阅读顺序",
        [
            "先看 dotplot：确认哪些通路进入 top10，以及通路名称是否与休眠、炎症、ECM remodeling、代谢重编程、免疫微环境或 TF 靶集相关。",
            "再看结果表：用 NES 判断富集方向，用 p.adjust/qvalue 判断统计证据，用 leading_edge 定位驱动富集的核心基因。",
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
            "每个 MSigDB 类别先展示 dotplot 图，再展示同一 analysis/geneset 对应的 GSEA 结果表。图中 top 通路可直接与下一页 NES、p.adjust、qvalue 和 leading_edge 统计量对应查看。",
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
                        "下一页是同一 GSEA 设计的结果表，请用 NES、p.adjust、qvalue 与 leading_edge 验证图中通路是否具有可靠统计和可追溯核心基因。",
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
                    "ID/Description 为通路编号与名称；setSize 为该通路基因集大小；enrichmentScore 和 NES 分别为原始富集分数与标准化富集分数。",
                    "本报告与 06 号运算统一采用 p.adjust < 0.05 作为 GSEA 图表展示阈值；qvalue 可作为额外稳健性参考。",
                    "rank 表示富集峰在 ranked gene list 中的位置；leading_edge 概括贡献最大的核心基因，是后续连接 DEG、TF 和药物重定位的关键字段。",
                    "正 NES 通常解释为 LRC 方向富集，负 NES 通常解释为 BULK 方向富集。",
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
            "Consensus_Rank 为综合排序；TF 为转录因子 symbol；Required_Methods/Source_Methods 显示该候选来自哪些方法组合。",
            "Source_Method_Count 表示支持该 TF 的方法数；Mean_Selected_Rank 与 Best_Selected_Rank 反映跨方法排序基础；CheA3_Library_Count 表示 ChEA3 多证据库支持数。",
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
        r"\newcommand{\ReportBodyFont}{\fontsize{5.65pt}{6.65pt}\selectfont}",
        r"\setbeamerfont{normal text}{size=\scriptsize}",
        r"\setbeamerfont{frametitle}{size=\large,series=\bfseries}",
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
        r"  \begin{frame}[plain]%",
        r"    \vfill%",
        r"    {\Large\bfseries #1}\par\vspace{0.45em}%",
        r"    {\normalsize #2}\par\vspace{0.8em}%",
        r"    \begin{minipage}{0.88\textwidth}\RaggedRight\scriptsize #3\end{minipage}%",
        r"    \vfill%",
        r"  \end{frame}%",
        r"}",
        "",
        r"\newcommand{\TextFrame}[2]{%",
        r"  \begin{frame}{#1}%",
        r"    \RaggedRight\ReportBodyFont #2%",
        r"  \end{frame}%",
        r"}",
        "",
        r"\newcommand{\ResultFigureFrame}[3]{%",
        r"  \begin{frame}[plain]%",
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
        r"  \begin{frame}[plain]%",
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
        r"  \begin{frame}[plain]{#1}%",
        r"    \vspace{-0.55em}%",
        r"    \begin{minipage}[t]{\textwidth}\RaggedRight\ReportBodyFont #4\end{minipage}%",
        r"    \par\vspace{0.08em}%",
        r"    \IfFileExists{\CRLMROOT/\detokenize{#3}}{%",
        r"      \begingroup\centering%",
        r"      \begin{adjustbox}{max width=0.998\textwidth,max totalheight=0.792\textheight,keepaspectratio}%",
        r"        \input{\CRLMROOT/\detokenize{#3}}%",
        r"      \end{adjustbox}%",
        r"      \par\endgroup%",
        r"    }{\MissingBox{#3}}%",
        r"  \end{frame}%",
        r"}",
        "",
        r"\newcommand{\ResultTwoTableFrame}[7]{%",
        r"  \begin{frame}[plain]{#1}%",
        r"    \vspace{-0.55em}%",
        r"    \noindent\begin{minipage}[t]{0.498\textwidth}\vspace{0pt}%",
        r"      \RaggedRight\ReportBodyFont #4%",
        r"      \par\vspace{0.08em}%",
        r"      \IfFileExists{\CRLMROOT/\detokenize{#3}}{%",
        r"        \begingroup\centering%",
        r"        \begin{adjustbox}{max width=0.995\linewidth,max totalheight=0.665\textheight,keepaspectratio}%",
        r"          \input{\CRLMROOT/\detokenize{#3}}%",
        r"        \end{adjustbox}%",
        r"        \par\endgroup%",
        r"      }{\MissingBox{#3}}%",
        r"    \end{minipage}%",
        r"    \hfill%",
        r"    \begin{minipage}[t]{0.498\textwidth}\vspace{0pt}%",
        r"      \RaggedRight\ReportBodyFont #7%",
        r"      \par\vspace{0.08em}%",
        r"      \IfFileExists{\CRLMROOT/\detokenize{#6}}{%",
        r"        \begingroup\centering%",
        r"        \begin{adjustbox}{max width=0.995\linewidth,max totalheight=0.665\textheight,keepaspectratio}%",
        r"          \input{\CRLMROOT/\detokenize{#6}}%",
        r"        \end{adjustbox}%",
        r"        \par\endgroup%",
        r"      }{\MissingBox{#6}}%",
        r"    \end{minipage}%",
        r"  \end{frame}%",
        r"}",
        "",
        r"\newcommand{\ResultStackedTwoTableFrame}[7]{%",
        r"  \begin{frame}[plain]{#1}%",
        r"    \vspace{-0.55em}%",
        r"    \begin{minipage}[t]{\textwidth}\vspace{0pt}%",
        r"      \RaggedRight\ReportBodyFont #4%",
        r"      \par\vspace{0.04em}%",
        r"      \IfFileExists{\CRLMROOT/\detokenize{#3}}{%",
        r"        \begingroup\centering%",
        r"        \begin{adjustbox}{max width=0.998\textwidth,max totalheight=0.305\textheight,keepaspectratio}%",
        r"          \input{\CRLMROOT/\detokenize{#3}}%",
        r"        \end{adjustbox}%",
        r"        \par\endgroup%",
        r"      }{\MissingBox{#3}}%",
        r"    \end{minipage}%",
        r"    \par\vspace{0.10em}%",
        r"    \begin{minipage}[t]{\textwidth}\vspace{0pt}%",
        r"      \RaggedRight\ReportBodyFont #7%",
        r"      \par\vspace{0.04em}%",
        r"      \IfFileExists{\CRLMROOT/\detokenize{#6}}{%",
        r"        \begingroup\centering%",
        r"        \begin{adjustbox}{max width=0.998\textwidth,max totalheight=0.305\textheight,keepaspectratio}%",
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
    section_names.append(build_01())
    section_names.append(build_02())
    section_names.append(build_04())
    section_names.append(build_05())
    section_names.extend(build_06())
    section_names.extend(build_07())
    section_names.extend(build_08())
    section_names.extend(build_09())
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
