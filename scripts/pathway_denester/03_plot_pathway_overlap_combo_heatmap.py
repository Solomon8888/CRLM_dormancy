#!/usr/bin/env python3
"""Draw PathwayDenester-style pathway-overlap combo heatmaps.

The WeChat/Paper example uses a two-panel Matplotlib layout:
top bars for -log10(p value), bottom heatmap for pairwise pathway overlap.
For GSEA outputs in this project, overlap is computed from each term's
core_enrichment / PathwayDenester intersection genes.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
except ImportError as exc:  # pragma: no cover - helpful runtime message
    raise SystemExit(
        "matplotlib is required for combo heatmaps. "
        "Install it in the PathwayDenester Python environment with: "
        "python -m pip install matplotlib"
    ) from exc


RESULTS_ROOT = Path("results")
SUMMARY_FILE = RESULTS_ROOT / "PathwayDenester_summary" / "summary.csv"
TOP_N = 20

TEXT_COLOR = "black"
TEXT_FONT_FAMILY = "Helvetica"
BAR_COLOR = "#334A9B"
BAR_EDGE_COLOR = "#1E2F69"
HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "pathway_overlap_blue",
    ["#FFFFFF", "#D8D7EC", "#8E8CC6", BAR_COLOR],
)


@dataclass(frozen=True)
class ComboHeatmapTask:
    data_type: str
    dataset_id: str
    plot_category: str
    analysis_name: str
    geneset_name: str
    input_file: Path
    result_file: Path
    plot_dir: Path
    table_dir: Path


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


def make_display_label(label: object, max_chars: int) -> str:
    label = str(label).replace("_", " ").strip()
    label = re.sub(r"\s+", " ", label)
    if len(label) <= max_chars:
        return label
    return label[: max(1, max_chars - 3)].rstrip() + "..."


def make_unique_labels(labels: Iterable[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for label in labels:
        count = seen.get(label, 0) + 1
        seen[label] = count
        out.append(label if count == 1 else f"{label} ({count})")
    return out


def compute_jaccard_matrix(gene_sets: list[set[str]]) -> np.ndarray:
    term_count = len(gene_sets)
    matrix = np.zeros((term_count, term_count), dtype=float)
    for i, genes_i in enumerate(gene_sets):
        for j, genes_j in enumerate(gene_sets):
            union_size = len(genes_i | genes_j)
            matrix[i, j] = 0.0 if union_size == 0 else len(genes_i & genes_j) / union_size
    return matrix


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

    input_dat = input_dat.copy()
    input_dat["p_value"] = pd.to_numeric(input_dat["p_value"], errors="coerce")
    if "intersection_size" not in input_dat.columns:
        input_dat["intersection_size"] = np.nan
    input_dat["Input_Order"] = np.arange(len(input_dat))
    input_dat["Genes"] = input_dat["intersection"].map(split_intersection_genes)
    input_dat = input_dat[
        input_dat["term_name"].notna()
        & input_dat["term_name"].astype(str).ne("")
        & np.isfinite(input_dat["p_value"])
        & input_dat["Genes"].map(len).gt(0)
    ].copy()
    input_dat = input_dat.sort_values(["p_value", "Input_Order"], kind="mergesort").head(top_n)

    if len(input_dat) < 2:
        return input_dat, np.zeros((0, 0), dtype=float)

    input_dat["Neg_Log10_P"] = safe_neglog10(input_dat["p_value"])
    input_dat["Pathway_Label"] = make_unique_labels(
        make_display_label(label, max_chars=56) for label in input_dat["term_name"]
    )

    gene_sets = [set(genes) for genes in input_dat["Genes"]]
    matrix = compute_jaccard_matrix(gene_sets)
    return input_dat.reset_index(drop=True), matrix


def get_figure_size(term_count: int) -> tuple[float, float]:
    width = min(max(9.6, term_count * 0.43 + 4.6), 14.5)
    height = min(max(8.2, term_count * 0.45 + 2.6), 13.8)
    return width, height


def apply_axis_text_style(ax: plt.Axes) -> None:
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily(TEXT_FONT_FAMILY)
        label.set_fontweight("bold")
        label.set_color(TEXT_COLOR)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)
        spine.set_linewidth(1.0)


def draw_empty_plot(task: ComboHeatmapTask, pdf_file: Path, png_file: Path, message: str) -> None:
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
    fig.savefig(pdf_file, bbox_inches="tight")
    fig.savefig(png_file, bbox_inches="tight", dpi=300)
    plt.close(fig)


def draw_combo_heatmap(
    task: ComboHeatmapTask,
    terms: pd.DataFrame,
    matrix: np.ndarray,
    pdf_file: Path,
    png_file: Path,
) -> None:
    term_count = len(terms)
    if term_count < 2 or matrix.size == 0:
        draw_empty_plot(
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
    width, height = get_figure_size(term_count)

    fig = plt.figure(figsize=(width, height))
    grid = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=[1.05, 5.5],
        left=0.17,
        right=0.74,
        bottom=0.20,
        top=0.92,
        hspace=0.06,
    )
    ax_bar = fig.add_subplot(grid[0, 0])
    ax_heat = fig.add_subplot(grid[1, 0])

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
    ax_bar.set_ylim(0, bar_top * 1.08 if bar_top > 0 else 1.0)
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
        fontsize=13,
        fontfamily=TEXT_FONT_FAMILY,
        fontweight="bold",
        color=TEXT_COLOR,
    )
    apply_axis_text_style(ax_bar)

    im = ax_heat.imshow(matrix, cmap=HEATMAP_CMAP, vmin=0, vmax=1, aspect="auto")
    ax_heat.set_xlim(-0.5, term_count - 0.5)
    ax_heat.set_ylim(term_count - 0.5, -0.5)
    ax_heat.set_xticks(x)
    ax_heat.set_xticklabels(
        [make_display_label(label, max_chars=42) for label in labels],
        rotation=-62,
        ha="left",
        va="top",
        rotation_mode="anchor",
        fontsize=8,
    )
    ax_heat.set_yticks(x)
    ax_heat.set_yticklabels(labels, fontsize=8)
    ax_heat.yaxis.tick_right()
    ax_heat.tick_params(axis="both", width=0, length=0, colors=TEXT_COLOR)

    ax_heat.set_xticks(np.arange(-0.5, term_count, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, term_count, 1), minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=0.7)
    ax_heat.tick_params(which="minor", bottom=False, left=False)
    apply_axis_text_style(ax_heat)

    heat_pos = ax_heat.get_position()
    cbar_width = 0.020
    cbar_height = min(0.18, heat_pos.height * 0.38)
    cbar_left = max(0.055, heat_pos.x0 - 0.070)
    cbar_bottom = heat_pos.y0 + heat_pos.height * 0.39
    cax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_ticks([0, 0.5, 1])
    cbar.set_ticklabels(["0", "0.5", "1"])
    cbar.ax.tick_params(labelsize=9, width=1.0, length=3, colors=TEXT_COLOR)
    cbar.outline.set_visible(False)
    for label in cbar.ax.get_yticklabels():
        label.set_fontfamily(TEXT_FONT_FAMILY)
        label.set_fontweight("bold")
        label.set_color(TEXT_COLOR)

    fig.text(
        cbar_left - 0.030,
        cbar_bottom + cbar_height / 2,
        "Pathway Overlap",
        ha="center",
        va="center",
        rotation=90,
        fontsize=13,
        fontfamily=TEXT_FONT_FAMILY,
        fontweight="bold",
        color=TEXT_COLOR,
    )

    fig.savefig(pdf_file, bbox_inches="tight")
    fig.savefig(png_file, bbox_inches="tight", dpi=300)
    plt.close(fig)


def write_tables(task: ComboHeatmapTask, terms: pd.DataFrame, matrix: np.ndarray) -> tuple[Path, Path]:
    task.table_dir.mkdir(parents=True, exist_ok=True)
    matrix_file = task.table_dir / "pathway_overlap_jaccard_matrix.csv"
    terms_file = task.table_dir / "pathway_overlap_heatmap_terms.csv"

    labels = terms["Pathway_Label"].tolist() if "Pathway_Label" in terms.columns else []
    if matrix.size == 0 or not labels:
        pd.DataFrame().to_csv(matrix_file, index=False)
    else:
        matrix_dat = pd.DataFrame(matrix, columns=labels)
        matrix_dat.insert(0, "Pathway", labels)
        matrix_dat.to_csv(matrix_file, index=False)

    out_terms = terms.drop(columns=["Genes"], errors="ignore").copy()
    out_terms.insert(0, "Heatmap_Rank", np.arange(1, len(out_terms) + 1))
    out_terms.to_csv(terms_file, index=False)
    return matrix_file, terms_file


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
            return {
                **base_record,
                "Terms_Plotted": int(len(terms)),
                "Matrix_File": str(task.table_dir / "pathway_overlap_jaccard_matrix.csv"),
                "Term_Metadata_File": str(task.table_dir / "pathway_overlap_heatmap_terms.csv"),
                "Status": "Skipped: existing plot",
                "Error_Message": "",
            }

        terms, matrix = prepare_combo_heatmap_data(task, top_n=top_n)
        matrix_file, terms_file = write_tables(task, terms, matrix)
        draw_combo_heatmap(
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
            "Status": "OK",
            "Error_Message": "",
        }
    except Exception as exc:  # keep the batch running, summarize failures
        return {
            **base_record,
            "Terms_Plotted": 0,
            "Matrix_File": "",
            "Term_Metadata_File": "",
            "Status": "ERROR",
            "Error_Message": str(exc),
        }


def write_summary_tables(summary: pd.DataFrame, results_root: Path) -> None:
    global_dir = results_root / "PathwayDenester_combo_heatmap_summary"
    global_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(global_dir / "summary.csv", index=False)

    for (data_type, dataset_id), dat in summary.groupby(["Data_Type", "Dataset_ID"], sort=True):
        out_dir = results_root / data_type / dataset_id / "tables" / "PathwayDenester_combo_heatmap_summary"
        out_dir.mkdir(parents=True, exist_ok=True)
        dat.to_csv(out_dir / "summary.csv", index=False)


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
        default=1,
        help="Accepted for command-line consistency; plotting runs sequentially.",
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

    records: list[dict[str, object]] = []
    for index, task in enumerate(tasks, start=1):
        records.append(plot_one_task(task, top_n=args.top_n, refresh=args.refresh))
        if index == 1 or index % max(1, min(20, len(tasks) // 10)) == 0 or index == len(tasks):
            print(f"Completed {index}/{len(tasks)}")

    summary = pd.DataFrame(records).sort_values(
        ["Data_Type", "Dataset_ID", "Plot_Category", "Analysis_Name", "GeneSet_Name"]
    )
    write_summary_tables(summary, results_root)

    ok_count = int(summary["Status"].astype(str).eq("OK").sum())
    skipped_count = int(summary["Status"].astype(str).str.startswith("Skipped").sum())
    error_count = int(summary["Status"].astype(str).eq("ERROR").sum())
    print("\nPathwayDenester combo heatmap plotting finished.")
    print(f"OK: {ok_count}")
    print(f"Skipped: {skipped_count}")
    print(f"Errors: {error_count}")
    print(f"Global summary: {results_root / 'PathwayDenester_combo_heatmap_summary' / 'summary.csv'}")
    return 1 if error_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
