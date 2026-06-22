#!/usr/bin/env python3
"""Draw PathwayDenester-style pathway-overlap combo heatmaps.

The WeChat/Paper example uses a two-panel Matplotlib layout:
top bars for -log10(p value), bottom heatmap for pairwise pathway overlap.
For GSEA outputs in this project, plotted terms are selected from the full
GSEA/PathwayDenester input before keep/exclude filtering, annotated with the
PathwayDenester "Filtered" status, and clustered by pathway overlap computed
from the full GMT gene sets.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import logging
import os
import re
import sys
import textwrap
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    from matplotlib.patches import Patch
except ImportError as exc:  # pragma: no cover - helpful runtime message
    raise SystemExit(
        "matplotlib is required for combo heatmaps. "
        "Install it in the PathwayDenester Python environment with: "
        "python -m pip install matplotlib"
    ) from exc


RESULTS_ROOT = Path("results")
SUMMARY_FILE = RESULTS_ROOT / "PathwayDenester_summary" / "summary.csv"
TOP_N = 40
PATHWAY_LABEL_WRAP_WIDTH = 45
PATHWAY_LABEL_PREFIXES = (
    "HALLMARK",
    "REACTOME",
    "KEGG_MEDICUS",
    "KEGG_LEGACY",
    "KEGG",
    "BIOCARTA",
    "WIKIPATHWAYS",
    "WP",
    "GO_BP",
    "GO_CC",
    "GO_MF",
    "GOBP",
    "GOCC",
    "GOMF",
    "HPO",
    "HP",
    "TFT_LEGACY",
    "TFT",
    "GTRD",
)
OVERLAP_METHOD = "row_pathway_coverage_intersection_over_row_pathway_size"
CLUSTERING_METHOD = "complete_linkage_on_one_minus_symmetrized_overlap"
TERM_SELECTION_METHOD = "top_n_full_gsea_input_by_p_value_before_pathwaydenester_filtering"

TEXT_COLOR = "black"
TEXT_FONT_FAMILY = "Helvetica"
TEXT_FONT_WEIGHT = "bold"
KEEP_COLOR = "#7D8BC5"
EXCLUDE_COLOR = "#F9C4DA"
UNKNOWN_COLOR = "#D8D8D8"
BAR_COLOR = "#33469D"
BAR_EDGE_COLOR = "#1E2F69"
HEATMAP_LOW_COLOR = "#FFFFFF"
HEATMAP_MID_LOW_COLOR = "#ECECF6"
HEATMAP_MID_COLOR = "#C8C8E0"
HEATMAP_MID_HIGH_COLOR = "#7D8BC5"
HEATMAP_HIGH_COLOR = "#33469D"
HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "pathway_overlap_blue",
    [
        HEATMAP_LOW_COLOR,
        HEATMAP_MID_LOW_COLOR,
        HEATMAP_MID_COLOR,
        HEATMAP_MID_HIGH_COLOR,
        HEATMAP_HIGH_COLOR,
    ],
)
FILTER_STATUS_COLORS = {
    "exclude": EXCLUDE_COLOR,
    "keep": KEEP_COLOR,
    "unknown": UNKNOWN_COLOR,
}

HEATMAP_CELL_SIZE_IN = 0.255
HEATMAP_MIN_SIDE_IN = 6.6
HEATMAP_MAX_SIDE_IN = 10.8
TOP_BAR_HEIGHT_RATIO = 0.18
TOP_BAR_GAP_IN = 0.12
STATUS_STRIP_WIDTH_IN = 0.32
STATUS_STRIP_GAP_IN = 0.11
COLORBAR_WIDTH_IN = 0.22
COLORBAR_HEIGHT_RATIO = 0.30
LEGEND_LEFT_IN = 0.20
LEFT_PANEL_WIDTH_IN = 2.25
RIGHT_LABEL_MIN_WIDTH_IN = 2.80
RIGHT_LABEL_MAX_WIDTH_IN = 5.80
FIGURE_RIGHT_MARGIN_IN = 0.40
FIGURE_TOP_MARGIN_IN = 0.28
FIGURE_BOTTOM_MIN_IN = 1.70
FIGURE_BOTTOM_LINE_HEIGHT_IN = 0.22
FIGURE_BOTTOM_CHAR_FACTOR_IN = 0.030
LABEL_FONT_SIZE_MIN = 5.9
LABEL_FONT_SIZE_MAX = 8.4
GRID_LINE_WIDTH = 0.48
AXIS_LINE_WIDTH = 1.0
PNG_DPI = 300


def configure_quiet_logging() -> None:
    warnings.filterwarnings("ignore")
    logging.getLogger("fontTools").setLevel(logging.ERROR)
    logging.getLogger("matplotlib").setLevel(logging.ERROR)


configure_quiet_logging()


@dataclass(frozen=True)
class ComboHeatmapTask:
    data_type: str
    dataset_id: str
    plot_category: str
    analysis_name: str
    geneset_name: str
    input_file: Path
    result_file: Path
    gmt_file: Path
    plot_dir: Path
    table_dir: Path


@dataclass(frozen=True)
class ComboHeatmapLayout:
    figure_width: float
    figure_height: float
    heatmap_left: float
    heatmap_bottom: float
    heatmap_side: float
    status_left: float
    status_width: float
    top_bar_bottom: float
    top_bar_height: float
    colorbar_left: float
    colorbar_bottom: float
    colorbar_width: float
    colorbar_height: float
    legend_left: float
    legend_bottom: float
    right_label_width: float
    bottom_margin: float
    label_font_size: float


def sanitize_file_name(value: object, default: str = "analysis") -> str:
    value = str(value).strip()
    if not value or value.lower() == "nan":
        value = default
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = re.sub(r"_+", "_", value)
    value = value.strip("_")
    return value or default


def parse_csv_arg(value: str) -> set[str] | None:
    value = (value or "").strip()
    if not value or value.lower() == "all":
        return None
    return {item.strip() for item in value.split(",") if item.strip()}


def split_intersection_genes(value: object) -> list[str]:
    if pd.isna(value):
        return []
    genes = re.split(r"[/,;]\s*", str(value).strip())
    seen: set[str] = set()
    out: list[str] = []
    for gene in genes:
        gene = gene.strip()
        if gene and gene not in seen:
            seen.add(gene)
            out.append(gene)
    return out


def safe_neglog10(p_values: pd.Series) -> np.ndarray:
    numeric = pd.to_numeric(p_values, errors="coerce").astype(float)
    positive = numeric[np.isfinite(numeric) & (numeric > 0)]
    floor = 1e-300 if positive.empty else max(float(positive.min()) * 0.1, 1e-300)
    safe = numeric.where(np.isfinite(numeric) & (numeric > 0), floor)
    safe = safe.clip(lower=1e-300)
    return -np.log10(safe.to_numpy(dtype=float))


def normalize_filter_status(value: object) -> str:
    status = str(value).strip().lower()
    if status in {"exclude", "keep"}:
        return status
    return "unknown"


def remove_pathway_prefix(label: object) -> str:
    label = str(label).strip()
    prefix_pattern = "|".join(re.escape(prefix) for prefix in PATHWAY_LABEL_PREFIXES)
    label = re.sub(rf"^({prefix_pattern})[_:\-\s]+", "", label, flags=re.IGNORECASE)
    label = re.sub(r"\s+", " ", label)
    label = label.replace("_", " ").strip()
    label = re.sub(rf"^({prefix_pattern})\s+", "", label, flags=re.IGNORECASE)
    return label or str(label)


def wrap_pathway_label(label: object, width: int = PATHWAY_LABEL_WRAP_WIDTH) -> str:
    label = remove_pathway_prefix(label)
    if len(label) <= width:
        return label

    wrapped = textwrap.wrap(
        label,
        width=width,
        break_long_words=False,
        break_on_hyphens=False,
    )
    if not wrapped or any(len(line) > width for line in wrapped):
        wrapped = textwrap.wrap(
            label,
            width=width,
            break_long_words=True,
            break_on_hyphens=False,
        )
    return "\n".join(wrapped)


def make_unique_labels(labels: Iterable[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for label in labels:
        count = seen.get(label, 0) + 1
        seen[label] = count
        out.append(label if count == 1 else f"{label} ({count})")
    return out


def compute_pathway_overlap_matrix(gene_sets: list[set[str]]) -> np.ndarray:
    term_count = len(gene_sets)
    matrix = np.zeros((term_count, term_count), dtype=float)
    for i, genes_i in enumerate(gene_sets):
        pathway_size = len(genes_i)
        for j, genes_j in enumerate(gene_sets):
            matrix[i, j] = 0.0 if pathway_size == 0 else len(genes_i & genes_j) / pathway_size
    return matrix


def get_cluster_distance(cluster_a: list[int], cluster_b: list[int], distance: np.ndarray) -> float:
    return max(float(distance[i, j]) for i in cluster_a for j in cluster_b)


def get_overlap_cluster_order(overlap_matrix: np.ndarray) -> list[int]:
    term_count = overlap_matrix.shape[0]
    if term_count < 3:
        return list(range(term_count))

    similarity = np.maximum(overlap_matrix, overlap_matrix.T)
    np.fill_diagonal(similarity, 1.0)
    distance = 1.0 - similarity
    np.fill_diagonal(distance, 0.0)

    clusters: list[list[int]] = [[i] for i in range(term_count)]
    while len(clusters) > 1:
        best_pair: tuple[int, int] | None = None
        best_distance = float("inf")
        for i in range(len(clusters) - 1):
            for j in range(i + 1, len(clusters)):
                candidate_distance = get_cluster_distance(clusters[i], clusters[j], distance)
                if candidate_distance < best_distance:
                    best_distance = candidate_distance
                    best_pair = (i, j)

        if best_pair is None:
            break

        i, j = best_pair
        clusters[i] = clusters[i] + clusters[j]
        del clusters[j]

    return clusters[0] if clusters else list(range(term_count))


def load_gmt_gene_sets(gmt_file: Path) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    if not gmt_file.exists():
        return gene_sets

    with gmt_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term_id = parts[0].strip()
            genes = {gene.strip() for gene in parts[2:] if gene.strip()}
            if term_id and genes:
                gene_sets[term_id] = genes
    return gene_sets


def load_filter_status(task: ComboHeatmapTask) -> pd.DataFrame:
    result_dat = pd.read_csv(task.result_file)
    if "Pathway ID" not in result_dat.columns:
        raise ValueError(f"{task.result_file} is missing Pathway ID column")

    optional_columns = [
        "Pathway ID",
        "Filtered",
        "Versus",
        "Versus Name",
        "Result",
        "Reciprocal pvalue",
        "DEGs in Intersection",
        "Is Reciprocal",
    ]
    optional_columns = [col for col in optional_columns if col in result_dat.columns]
    status_dat = result_dat[optional_columns].copy()
    status_dat["term_id"] = status_dat["Pathway ID"].astype(str)
    status_dat["Filtered"] = status_dat.get("Filtered", "unknown")
    status_dat["Filtered"] = status_dat["Filtered"].map(normalize_filter_status)
    status_dat = status_dat.drop_duplicates("term_id", keep="first")
    return status_dat.drop(columns=["Pathway ID"], errors="ignore")


def discover_tasks(
    summary_file: Path,
    results_root: Path,
    datasets: set[str] | None,
    genesets: set[str] | None,
    categories: set[str] | None,
) -> list[ComboHeatmapTask]:
    if not summary_file.exists():
        raise FileNotFoundError(f"Missing PathwayDenester summary: {summary_file}")

    summary = pd.read_csv(summary_file)
    keep_status = summary["Status"].astype(str).isin(["OK", "Skipped: existing result"])
    summary = summary.loc[keep_status].copy()

    if datasets is not None:
        summary = summary[summary["Dataset_ID"].astype(str).isin(datasets)]
    if genesets is not None:
        summary = summary[summary["GeneSet_Name"].astype(str).isin(genesets)]
    if categories is not None:
        summary = summary[summary["Plot_Category"].astype(str).isin(categories)]

    tasks: list[ComboHeatmapTask] = []
    for _, row in summary.sort_values(
        ["Data_Type", "Dataset_ID", "Plot_Category", "Analysis_Name", "GeneSet_Name"]
    ).iterrows():
        data_type = str(row["Data_Type"])
        dataset_id = str(row["Dataset_ID"])
        plot_category = str(row["Plot_Category"])
        analysis_name = str(row["Analysis_Name"])
        geneset_name = str(row["GeneSet_Name"])

        table_dir = Path(str(row["PathwayDenester_Result_File"])).parent
        plot_dir = (
            results_root
            / data_type
            / dataset_id
            / "plots"
            / "PathwayDenester"
            / sanitize_file_name(plot_category)
            / sanitize_file_name(analysis_name)
            / sanitize_file_name(geneset_name)
        )
        tasks.append(
            ComboHeatmapTask(
                data_type=data_type,
                dataset_id=dataset_id,
                plot_category=plot_category,
                analysis_name=analysis_name,
                geneset_name=geneset_name,
                input_file=Path(str(row["PathwayDenester_Input_File"])),
                result_file=Path(str(row["PathwayDenester_Result_File"])),
                gmt_file=Path(str(row["GMT_File"])),
                plot_dir=plot_dir,
                table_dir=table_dir,
            )
        )

    return tasks


def prepare_combo_heatmap_data(task: ComboHeatmapTask, top_n: int) -> tuple[pd.DataFrame, np.ndarray]:
    input_dat = pd.read_csv(task.input_file, sep="\t")
    required = {"term_id", "term_name", "p_value", "intersection"}
    missing = sorted(required - set(input_dat.columns))
    if missing:
        raise ValueError(f"{task.input_file} is missing columns: {', '.join(missing)}")

    gmt_gene_sets = load_gmt_gene_sets(task.gmt_file)
    input_dat = input_dat.copy()
    input_dat["term_id"] = input_dat["term_id"].astype(str)
    input_dat["p_value"] = pd.to_numeric(input_dat["p_value"], errors="coerce")
    if "intersection_size" not in input_dat.columns:
        input_dat["intersection_size"] = np.nan
    input_dat["Input_Order"] = np.arange(len(input_dat))
    input_dat["Selected_Genes"] = input_dat["intersection"].map(split_intersection_genes)
    input_dat["Pathway_Genes"] = input_dat.apply(
        lambda row: sorted(gmt_gene_sets.get(row["term_id"], set(row["Selected_Genes"]))),
        axis=1,
    )
    input_dat["Pathway_Gene_Source"] = input_dat["term_id"].map(
        lambda term_id: "GMT" if term_id in gmt_gene_sets else "selected_genes_fallback"
    )
    input_dat["Pathway_Gene_Count"] = input_dat["Pathway_Genes"].map(len)

    status_dat = load_filter_status(task)
    input_dat = input_dat.merge(status_dat, on="term_id", how="left")
    input_dat["Filtered"] = input_dat["Filtered"].map(normalize_filter_status)

    input_dat = input_dat[
        input_dat["term_name"].notna()
        & input_dat["term_name"].astype(str).ne("")
        & np.isfinite(input_dat["p_value"])
        & input_dat["Pathway_Genes"].map(len).gt(0)
    ].copy()
    input_dat = input_dat.sort_values(["p_value", "Input_Order"], kind="mergesort")
    input_dat["GSEA_Rank"] = np.arange(1, len(input_dat) + 1)
    if top_n > 0:
        input_dat = input_dat.head(top_n)

    if len(input_dat) < 2:
        return input_dat, np.zeros((0, 0), dtype=float)

    input_dat["Neg_Log10_P"] = safe_neglog10(input_dat["p_value"])
    input_dat["Pathway_Label_Clean"] = input_dat["term_name"].map(remove_pathway_prefix)
    input_dat["Pathway_Label"] = make_unique_labels(
        wrap_pathway_label(label, width=PATHWAY_LABEL_WRAP_WIDTH)
        for label in input_dat["Pathway_Label_Clean"]
    )

    gene_sets = [set(genes) for genes in input_dat["Pathway_Genes"]]
    matrix = compute_pathway_overlap_matrix(gene_sets)
    order = get_overlap_cluster_order(matrix)
    input_dat = input_dat.iloc[order, :].reset_index(drop=True)
    matrix = matrix[np.ix_(order, order)]
    input_dat["Heatmap_Cluster_Rank"] = np.arange(1, len(input_dat) + 1)
    return input_dat.reset_index(drop=True), matrix


def get_label_stats(labels: list[str]) -> tuple[int, int]:
    max_lines = 1
    max_line_length = 1
    for label in labels:
        lines = str(label).split("\n")
        max_lines = max(max_lines, len(lines))
        max_line_length = max(max_line_length, *(len(line) for line in lines))
    return max_lines, max_line_length


def get_label_font_size(term_count: int) -> float:
    if term_count <= 20:
        return LABEL_FONT_SIZE_MAX
    if term_count >= 50:
        return LABEL_FONT_SIZE_MIN
    fraction = (term_count - 20) / 30
    return LABEL_FONT_SIZE_MAX - fraction * (LABEL_FONT_SIZE_MAX - LABEL_FONT_SIZE_MIN)


def get_combo_heatmap_layout(term_count: int, labels: list[str]) -> ComboHeatmapLayout:
    label_max_lines, label_max_line_length = get_label_stats(labels)
    heatmap_side = min(
        max(term_count * HEATMAP_CELL_SIZE_IN, HEATMAP_MIN_SIDE_IN),
        HEATMAP_MAX_SIDE_IN,
    )
    top_bar_height = heatmap_side * TOP_BAR_HEIGHT_RATIO
    bottom_margin = max(
        FIGURE_BOTTOM_MIN_IN,
        0.48
        + label_max_lines * FIGURE_BOTTOM_LINE_HEIGHT_IN
        + min(label_max_line_length, PATHWAY_LABEL_WRAP_WIDTH) * FIGURE_BOTTOM_CHAR_FACTOR_IN,
    )
    right_label_width = min(
        max(
            RIGHT_LABEL_MIN_WIDTH_IN,
            0.65
            + min(label_max_line_length, PATHWAY_LABEL_WRAP_WIDTH) * 0.070
            + label_max_lines * 0.12,
        ),
        RIGHT_LABEL_MAX_WIDTH_IN,
    )

    status_left = LEFT_PANEL_WIDTH_IN
    heatmap_left = status_left + STATUS_STRIP_WIDTH_IN + STATUS_STRIP_GAP_IN
    heatmap_bottom = bottom_margin
    top_bar_bottom = heatmap_bottom + heatmap_side + TOP_BAR_GAP_IN
    colorbar_height = heatmap_side * COLORBAR_HEIGHT_RATIO
    colorbar_left = max(LEGEND_LEFT_IN + 0.55, status_left - 1.06)
    colorbar_bottom = heatmap_bottom + heatmap_side * 0.50
    figure_width = heatmap_left + heatmap_side + right_label_width + FIGURE_RIGHT_MARGIN_IN
    figure_height = bottom_margin + heatmap_side + TOP_BAR_GAP_IN + top_bar_height + FIGURE_TOP_MARGIN_IN

    return ComboHeatmapLayout(
        figure_width=figure_width,
        figure_height=figure_height,
        heatmap_left=heatmap_left,
        heatmap_bottom=heatmap_bottom,
        heatmap_side=heatmap_side,
        status_left=status_left,
        status_width=STATUS_STRIP_WIDTH_IN,
        top_bar_bottom=top_bar_bottom,
        top_bar_height=top_bar_height,
        colorbar_left=colorbar_left,
        colorbar_bottom=colorbar_bottom,
        colorbar_width=COLORBAR_WIDTH_IN,
        colorbar_height=colorbar_height,
        legend_left=LEGEND_LEFT_IN,
        legend_bottom=heatmap_bottom + 0.10,
        right_label_width=right_label_width,
        bottom_margin=bottom_margin,
        label_font_size=get_label_font_size(term_count),
    )


def normalize_axes_box(layout: ComboHeatmapLayout, left: float, bottom: float, width: float, height: float) -> list[float]:
    return [
        left / layout.figure_width,
        bottom / layout.figure_height,
        width / layout.figure_width,
        height / layout.figure_height,
    ]


def apply_axis_text_style(ax: plt.Axes) -> None:
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily(TEXT_FONT_FAMILY)
        label.set_fontweight(TEXT_FONT_WEIGHT)
        label.set_color(TEXT_COLOR)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)
        spine.set_linewidth(AXIS_LINE_WIDTH)


def save_figure_quietly(fig: plt.Figure, pdf_file: Path, png_file: Path) -> None:
    with open(os.devnull, "w", encoding="utf-8") as devnull:
        with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
            fig.savefig(pdf_file, bbox_inches="tight")
            fig.savefig(png_file, bbox_inches="tight", dpi=300)


def get_layout_parameter_record(layout: ComboHeatmapLayout | None = None) -> dict[str, object]:
    if layout is None:
        return {
            "Figure_Width_In": np.nan,
            "Figure_Height_In": np.nan,
            "Heatmap_Side_In": np.nan,
            "Top_Bar_Height_In": np.nan,
            "Status_Strip_Width_In": STATUS_STRIP_WIDTH_IN,
            "Bottom_Margin_In": np.nan,
            "Right_Label_Width_In": np.nan,
            "Label_Font_Size": np.nan,
        }

    return {
        "Figure_Width_In": round(layout.figure_width, 4),
        "Figure_Height_In": round(layout.figure_height, 4),
        "Heatmap_Side_In": round(layout.heatmap_side, 4),
        "Top_Bar_Height_In": round(layout.top_bar_height, 4),
        "Status_Strip_Width_In": round(layout.status_width, 4),
        "Bottom_Margin_In": round(layout.bottom_margin, 4),
        "Right_Label_Width_In": round(layout.right_label_width, 4),
        "Label_Font_Size": round(layout.label_font_size, 4),
    }


def draw_empty_plot(task: ComboHeatmapTask, pdf_file: Path, png_file: Path, message: str) -> dict[str, object]:
    pdf_file.parent.mkdir(parents=True, exist_ok=True)
    png_file.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        fontsize=13,
        fontfamily=TEXT_FONT_FAMILY,
        fontweight="bold",
        color=TEXT_COLOR,
    )
    ax.set_axis_off()
    save_figure_quietly(fig, pdf_file=pdf_file, png_file=png_file)
    plt.close(fig)
    return get_layout_parameter_record()


def draw_combo_heatmap(
    task: ComboHeatmapTask,
    terms: pd.DataFrame,
    matrix: np.ndarray,
    pdf_file: Path,
    png_file: Path,
) -> dict[str, object]:
    term_count = len(terms)
    if term_count < 2 or matrix.size == 0:
        return draw_empty_plot(
            task,
            pdf_file,
            png_file,
            "Fewer than two pathways with leading-edge genes",
        )
        return

    pdf_file.parent.mkdir(parents=True, exist_ok=True)
    png_file.parent.mkdir(parents=True, exist_ok=True)
    plt.rcParams["font.family"] = TEXT_FONT_FAMILY
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    labels = terms["Pathway_Label"].tolist()
    neglog = terms["Neg_Log10_P"].to_numpy(dtype=float)
    statuses = terms["Filtered"].map(normalize_filter_status).tolist()
    status_colors = [FILTER_STATUS_COLORS.get(status, UNKNOWN_COLOR) for status in statuses]
    layout = get_combo_heatmap_layout(term_count=term_count, labels=labels)

    fig = plt.figure(figsize=(layout.figure_width, layout.figure_height))
    ax_bar = fig.add_axes(
        normalize_axes_box(
            layout,
            layout.heatmap_left,
            layout.top_bar_bottom,
            layout.heatmap_side,
            layout.top_bar_height,
        )
    )
    ax_status = fig.add_axes(
        normalize_axes_box(
            layout,
            layout.status_left,
            layout.heatmap_bottom,
            layout.status_width,
            layout.heatmap_side,
        )
    )
    ax_heat = fig.add_axes(
        normalize_axes_box(
            layout,
            layout.heatmap_left,
            layout.heatmap_bottom,
            layout.heatmap_side,
            layout.heatmap_side,
        )
    )

    x = np.arange(term_count)
    ax_bar.bar(
        x,
        neglog,
        width=0.92,
        color=BAR_COLOR,
        edgecolor=BAR_EDGE_COLOR,
        linewidth=0.2,
    )
    bar_top = float(np.nanmax(neglog)) if np.isfinite(neglog).any() else 1.0
    ax_bar.set_ylim(0, bar_top * 1.20 if bar_top > 0 else 1.0)
    ax_bar.set_xlim(-0.5, term_count - 0.5)
    ax_bar.set_xticks([])
    ax_bar.tick_params(axis="y", labelsize=9, width=1.0, length=4, colors=TEXT_COLOR)
    ax_bar.text(
        1.01,
        0.55,
        "-log10(p-value)",
        transform=ax_bar.transAxes,
        ha="left",
        va="center",
        fontsize=13.2,
        fontfamily=TEXT_FONT_FAMILY,
        fontweight=TEXT_FONT_WEIGHT,
        color=TEXT_COLOR,
    )
    apply_axis_text_style(ax_bar)

    status_rgb = np.array([[to_rgb(color)] for color in status_colors])
    ax_status.imshow(status_rgb, aspect="auto")
    ax_status.set_ylim(term_count - 0.5, -0.5)
    ax_status.set_xticks([])
    ax_status.set_yticks([])
    for spine in ax_status.spines.values():
        spine.set_color(TEXT_COLOR)
        spine.set_linewidth(0.8)

    im = ax_heat.imshow(matrix, cmap=HEATMAP_CMAP, vmin=0, vmax=1, aspect="equal")
    ax_heat.set_xlim(-0.5, term_count - 0.5)
    ax_heat.set_ylim(term_count - 0.5, -0.5)
    ax_heat.set_xticks(x)
    ax_heat.set_xticklabels(
        labels,
        rotation=-62,
        ha="left",
        va="top",
        rotation_mode="anchor",
        fontsize=layout.label_font_size,
    )
    ax_heat.set_yticks(x)
    ax_heat.set_yticklabels(labels, fontsize=layout.label_font_size)
    ax_heat.yaxis.tick_right()
    ax_heat.tick_params(axis="both", width=0, length=0, colors=TEXT_COLOR, pad=3)

    ax_heat.set_xticks(np.arange(-0.5, term_count, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, term_count, 1), minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=GRID_LINE_WIDTH)
    ax_heat.tick_params(which="minor", bottom=False, left=False)
    apply_axis_text_style(ax_heat)
    for label in ax_heat.get_xticklabels() + ax_heat.get_yticklabels():
        label.set_linespacing(0.92)

    legend_handles = [
        Patch(facecolor=EXCLUDE_COLOR, edgecolor=TEXT_COLOR, label="exclude"),
        Patch(facecolor=KEEP_COLOR, edgecolor=TEXT_COLOR, label="keep"),
    ]
    if "unknown" in statuses:
        legend_handles.append(
            Patch(facecolor=UNKNOWN_COLOR, edgecolor=TEXT_COLOR, label="unknown")
        )
    legend = fig.legend(
        handles=legend_handles,
        title="filtered",
        loc="lower left",
        bbox_to_anchor=(
            layout.legend_left / layout.figure_width,
            layout.legend_bottom / layout.figure_height,
        ),
        frameon=False,
        fontsize=11,
        title_fontsize=12,
    )
    for text in legend.get_texts() + [legend.get_title()]:
        text.set_fontfamily(TEXT_FONT_FAMILY)
        text.set_fontweight(TEXT_FONT_WEIGHT)
        text.set_color(TEXT_COLOR)

    cax = fig.add_axes(
        normalize_axes_box(
            layout,
            layout.colorbar_left,
            layout.colorbar_bottom,
            layout.colorbar_width,
            layout.colorbar_height,
        )
    )
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_ticks([0, 0.5, 1])
    cbar.set_ticklabels(["0", "0.5", "1"])
    cbar.ax.tick_params(labelsize=9, width=1.0, length=3, colors=TEXT_COLOR)
    cbar.outline.set_visible(False)
    for label in cbar.ax.get_yticklabels():
        label.set_fontfamily(TEXT_FONT_FAMILY)
        label.set_fontweight(TEXT_FONT_WEIGHT)
        label.set_color(TEXT_COLOR)

    fig.text(
        (layout.colorbar_left - 0.37) / layout.figure_width,
        (layout.colorbar_bottom + layout.colorbar_height / 2) / layout.figure_height,
        "Pathway Overlap",
        ha="center",
        va="center",
        rotation=90,
        fontsize=13,
        fontfamily=TEXT_FONT_FAMILY,
        fontweight=TEXT_FONT_WEIGHT,
        color=TEXT_COLOR,
    )

    save_figure_quietly(fig, pdf_file=pdf_file, png_file=png_file)
    plt.close(fig)
    return get_layout_parameter_record(layout)


def write_tables(task: ComboHeatmapTask, terms: pd.DataFrame, matrix: np.ndarray) -> tuple[Path, Path]:
    task.table_dir.mkdir(parents=True, exist_ok=True)
    matrix_file = task.table_dir / "pathway_overlap_matrix.csv"
    terms_file = task.table_dir / "pathway_overlap_heatmap_terms.csv"
    legacy_jaccard_file = task.table_dir / "pathway_overlap_jaccard_matrix.csv"
    if legacy_jaccard_file.exists():
        legacy_jaccard_file.unlink()

    labels = terms["Pathway_Label"].tolist() if "Pathway_Label" in terms.columns else []
    if matrix.size == 0 or not labels:
        pd.DataFrame().to_csv(matrix_file, index=False)
    else:
        matrix_dat = pd.DataFrame(matrix, columns=labels)
        matrix_dat.insert(0, "Pathway", labels)
        matrix_dat.to_csv(matrix_file, index=False)

    out_terms = terms.drop(columns=["Selected_Genes", "Pathway_Genes"], errors="ignore").copy()
    out_terms.insert(0, "Heatmap_Rank", np.arange(1, len(out_terms) + 1))
    out_terms.to_csv(terms_file, index=False)
    return matrix_file, terms_file


def get_existing_layout_parameter_record(terms: pd.DataFrame) -> dict[str, object]:
    if "Pathway_Label" not in terms.columns or len(terms) < 2:
        return get_layout_parameter_record()
    labels = terms["Pathway_Label"].astype(str).tolist()
    layout = get_combo_heatmap_layout(term_count=len(terms), labels=labels)
    return get_layout_parameter_record(layout)


def plot_one_task(task: ComboHeatmapTask, top_n: int, refresh: bool) -> dict[str, object]:
    pdf_file = task.plot_dir / "pdf" / "pathway_overlap_heatmap.pdf"
    png_file = task.plot_dir / "png" / "pathway_overlap_heatmap.png"
    base_record = {
        "Data_Type": task.data_type,
        "Dataset_ID": task.dataset_id,
        "Plot_Category": task.plot_category,
        "Analysis_Name": task.analysis_name,
        "GeneSet_Name": task.geneset_name,
        "Plot_Type": "pathway_overlap_combo_heatmap",
        "PDF_File": str(pdf_file),
        "PNG_File": str(png_file),
    }

    try:
        if not task.input_file.exists():
            raise FileNotFoundError(f"Missing PathwayDenester input: {task.input_file}")
        if not task.result_file.exists():
            raise FileNotFoundError(f"Missing PathwayDenester result: {task.result_file}")

        if pdf_file.exists() and png_file.exists() and not refresh:
            terms = pd.read_csv(task.table_dir / "pathway_overlap_heatmap_terms.csv")
            plot_parameters = get_existing_layout_parameter_record(terms)
            return {
                **base_record,
                "Terms_Plotted": int(len(terms)),
                "Matrix_File": str(task.table_dir / "pathway_overlap_matrix.csv"),
                "Term_Metadata_File": str(task.table_dir / "pathway_overlap_heatmap_terms.csv"),
                **plot_parameters,
                "Status": "Skipped: existing plot",
                "Error_Message": "",
            }

        terms, matrix = prepare_combo_heatmap_data(task, top_n=top_n)
        matrix_file, terms_file = write_tables(task, terms, matrix)
        plot_parameters = draw_combo_heatmap(
            task=task,
            terms=terms,
            matrix=matrix,
            pdf_file=pdf_file,
            png_file=png_file,
        )

        return {
            **base_record,
            "Terms_Plotted": int(len(terms)),
            "Matrix_File": str(matrix_file),
            "Term_Metadata_File": str(terms_file),
            **plot_parameters,
            "Status": "OK",
            "Error_Message": "",
        }
    except Exception as exc:  # keep the batch running, summarize failures
        return {
            **base_record,
            "Terms_Plotted": 0,
            "Matrix_File": "",
            "Term_Metadata_File": "",
            **get_layout_parameter_record(),
            "Status": "ERROR",
            "Error_Message": str(exc),
        }


def plot_one_task_worker(payload: tuple[ComboHeatmapTask, int, bool]) -> dict[str, object]:
    configure_quiet_logging()
    task, top_n, refresh = payload
    return plot_one_task(task=task, top_n=top_n, refresh=refresh)


def get_method_parameter_table(top_n: int, workers: int) -> pd.DataFrame:
    parameters = {
        "top_n": top_n,
        "term_selection_method": TERM_SELECTION_METHOD,
        "overlap_method": OVERLAP_METHOD,
        "clustering_method": CLUSTERING_METHOD,
        "pathway_gene_source_priority": "GMT_full_gene_set_then_selected_genes_fallback",
        "p_value_transform": "-log10(raw_GSEA_pvalue)",
        "filter_annotation_source": "PathwayDenester_Filtered_column",
        "label_prefixes_removed": ",".join(PATHWAY_LABEL_PREFIXES),
        "label_wrap_width": PATHWAY_LABEL_WRAP_WIDTH,
        "text_font_family": TEXT_FONT_FAMILY,
        "text_font_weight": TEXT_FONT_WEIGHT,
        "bar_color": BAR_COLOR,
        "bar_edge_color": BAR_EDGE_COLOR,
        "heatmap_low_color": HEATMAP_LOW_COLOR,
        "heatmap_mid_low_color": HEATMAP_MID_LOW_COLOR,
        "heatmap_mid_color": HEATMAP_MID_COLOR,
        "heatmap_mid_high_color": HEATMAP_MID_HIGH_COLOR,
        "heatmap_high_color": HEATMAP_HIGH_COLOR,
        "keep_color": KEEP_COLOR,
        "exclude_color": EXCLUDE_COLOR,
        "unknown_filter_color": UNKNOWN_COLOR,
        "heatmap_cell_size_in": HEATMAP_CELL_SIZE_IN,
        "heatmap_min_side_in": HEATMAP_MIN_SIDE_IN,
        "heatmap_max_side_in": HEATMAP_MAX_SIDE_IN,
        "top_bar_height_ratio": TOP_BAR_HEIGHT_RATIO,
        "top_bar_gap_in": TOP_BAR_GAP_IN,
        "status_strip_width_in": STATUS_STRIP_WIDTH_IN,
        "status_strip_gap_in": STATUS_STRIP_GAP_IN,
        "left_panel_width_in": LEFT_PANEL_WIDTH_IN,
        "right_label_min_width_in": RIGHT_LABEL_MIN_WIDTH_IN,
        "right_label_max_width_in": RIGHT_LABEL_MAX_WIDTH_IN,
        "png_dpi": PNG_DPI,
        "workers": workers,
    }
    return pd.DataFrame(
        [
            {
                "Parameter": key,
                "Value": value,
            }
            for key, value in parameters.items()
        ]
    )


def add_common_summary_parameters(summary: pd.DataFrame, top_n: int, workers: int) -> pd.DataFrame:
    summary = summary.copy()
    summary["Top_N"] = top_n
    summary["Term_Selection_Method"] = TERM_SELECTION_METHOD
    summary["Overlap_Method"] = OVERLAP_METHOD
    summary["Clustering_Method"] = CLUSTERING_METHOD
    summary["PValue_Transform"] = "-log10(raw_GSEA_pvalue)"
    summary["Pathway_Gene_Source_Priority"] = "GMT_full_gene_set_then_selected_genes_fallback"
    summary["Label_Wrap_Width"] = PATHWAY_LABEL_WRAP_WIDTH
    summary["Text_Font_Family"] = TEXT_FONT_FAMILY
    summary["Text_Font_Weight"] = TEXT_FONT_WEIGHT
    summary["Keep_Color"] = KEEP_COLOR
    summary["Exclude_Color"] = EXCLUDE_COLOR
    summary["Heatmap_High_Color"] = HEATMAP_HIGH_COLOR
    summary["Bar_Color"] = BAR_COLOR
    summary["Workers_Used"] = workers
    return summary


def write_summary_tables(summary: pd.DataFrame, results_root: Path, top_n: int, workers: int) -> None:
    method_parameters = get_method_parameter_table(top_n=top_n, workers=workers)
    global_dir = results_root / "PathwayDenester_combo_heatmap_summary"
    global_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(global_dir / "summary.csv", index=False)
    method_parameters.to_csv(global_dir / "method_parameters.csv", index=False)

    for (data_type, dataset_id), dat in summary.groupby(["Data_Type", "Dataset_ID"], sort=True):
        out_dir = results_root / data_type / dataset_id / "tables" / "PathwayDenester_combo_heatmap_summary"
        out_dir.mkdir(parents=True, exist_ok=True)
        dat.to_csv(out_dir / "summary.csv", index=False)
        method_parameters.to_csv(out_dir / "method_parameters.csv", index=False)


def format_duration(seconds: float) -> str:
    seconds = max(0, int(seconds))
    minutes, sec = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return f"{hours:d}h{minutes:02d}m{sec:02d}s"
    if minutes:
        return f"{minutes:d}m{sec:02d}s"
    return f"{sec:d}s"


def build_progress_line(
    completed: int,
    total: int,
    ok_count: int,
    skipped_count: int,
    error_count: int,
    workers: int,
    start_time: float,
    width: int = 32,
) -> str:
    fraction = completed / total if total else 1.0
    filled = int(round(width * fraction))
    bar = "=" * filled + "." * (width - filled)
    elapsed = time.time() - start_time
    rate = completed / elapsed if elapsed > 0 else 0
    eta = (total - completed) / rate if rate > 0 else 0
    return (
        f"\rPathwayDenester combo heatmaps [{bar}] "
        f"{completed}/{total} {fraction * 100:5.1f}% | "
        f"OK {ok_count} Skip {skipped_count} Err {error_count} | "
        f"workers {workers} | elapsed {format_duration(elapsed)} | ETA {format_duration(eta)}"
    )


def run_plot_tasks(
    tasks: list[ComboHeatmapTask],
    top_n: int,
    refresh: bool,
    workers: int,
) -> pd.DataFrame:
    total = len(tasks)
    workers = max(1, min(workers, total))
    start_time = time.time()
    records: list[dict[str, object]] = []
    ok_count = 0
    skipped_count = 0
    error_count = 0
    last_progress_time = 0.0

    def emit_progress(force: bool = False) -> None:
        nonlocal last_progress_time
        now = time.time()
        if not force and now - last_progress_time < 1.0:
            return
        last_progress_time = now
        print(
            build_progress_line(
                completed=len(records),
                total=total,
                ok_count=ok_count,
                skipped_count=skipped_count,
                error_count=error_count,
                workers=workers,
                start_time=start_time,
            ),
            end="",
            flush=True,
        )

    print(
        build_progress_line(
            completed=0,
            total=total,
            ok_count=0,
            skipped_count=0,
            error_count=0,
            workers=workers,
            start_time=start_time,
        ),
        end="",
        flush=True,
    )
    last_progress_time = time.time()

    payloads = [(task, top_n, refresh) for task in tasks]
    if workers == 1:
        for payload in payloads:
            record = plot_one_task_worker(payload)
            records.append(record)
            status = str(record.get("Status", ""))
            ok_count += status == "OK"
            skipped_count += status.startswith("Skipped")
            error_count += status == "ERROR"
            emit_progress(force=len(records) == total)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
            future_to_task = {
                executor.submit(plot_one_task_worker, payload): payload[0]
                for payload in payloads
            }
            for future in concurrent.futures.as_completed(future_to_task):
                record = future.result()
                records.append(record)
                status = str(record.get("Status", ""))
                ok_count += status == "OK"
                skipped_count += status.startswith("Skipped")
                error_count += status == "ERROR"
                emit_progress(force=len(records) == total)

    print()
    return pd.DataFrame(records)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch-plot PathwayDenester pathway-overlap combo heatmaps."
    )
    parser.add_argument("--results-root", default=str(RESULTS_ROOT))
    parser.add_argument("--summary-file", default=str(SUMMARY_FILE))
    parser.add_argument("--datasets", default="all", help="Comma-separated dataset IDs or all.")
    parser.add_argument("--genesets", default="all", help="Comma-separated geneset output names or all.")
    parser.add_argument("--categories", default="all", help="Comma-separated plot categories or all.")
    parser.add_argument("--top-n", type=int, default=TOP_N)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument(
        "--clone-if-missing",
        action="store_true",
        help="Accepted for compatibility with 99_run_all_pathway_denester.R; ignored here.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=os.cpu_count() or 1,
        help="Parallel plotting workers. Use 0 or a negative value to use all logical CPU cores.",
    )
    return parser.parse_args(list(argv))


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    results_root = Path(args.results_root)
    tasks = discover_tasks(
        summary_file=Path(args.summary_file),
        results_root=results_root,
        datasets=parse_csv_arg(args.datasets),
        genesets=parse_csv_arg(args.genesets),
        categories=parse_csv_arg(args.categories),
    )
    if not tasks:
        raise SystemExit("No PathwayDenester rows were selected for combo heatmap plotting.")

    print("\nRunning PathwayDenester combo heatmap plotting...")
    print(f"Selected PathwayDenester results: {len(tasks)}")
    print(f"Top pathways per heatmap: {args.top_n}")
    workers = int(args.workers)
    if workers <= 0:
        workers = os.cpu_count() or 1
    workers = max(1, min(workers, len(tasks)))
    print(f"Workers: {workers}")

    summary = run_plot_tasks(
        tasks=tasks,
        top_n=args.top_n,
        refresh=args.refresh,
        workers=workers,
    ).sort_values(
        ["Data_Type", "Dataset_ID", "Plot_Category", "Analysis_Name", "GeneSet_Name"]
    )
    summary = add_common_summary_parameters(summary, top_n=args.top_n, workers=workers)
    write_summary_tables(summary, results_root, top_n=args.top_n, workers=workers)

    ok_count = int(summary["Status"].astype(str).eq("OK").sum())
    skipped_count = int(summary["Status"].astype(str).str.startswith("Skipped").sum())
    error_count = int(summary["Status"].astype(str).eq("ERROR").sum())
    print("\nPathwayDenester combo heatmap plotting finished.")
    print(f"OK: {ok_count}")
    print(f"Skipped: {skipped_count}")
    print(f"Errors: {error_count}")
    print(f"Global summary: {results_root / 'PathwayDenester_combo_heatmap_summary' / 'summary.csv'}")
    print(f"Method parameters: {results_root / 'PathwayDenester_combo_heatmap_summary' / 'method_parameters.csv'}")
    return 1 if error_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
