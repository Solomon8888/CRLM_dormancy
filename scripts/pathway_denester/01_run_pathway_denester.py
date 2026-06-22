#!/usr/bin/env python3
"""Run official PathwayDenester over all project GSEA result tables.

The official tool expects pathway enrichment rows with:
term_id, term_name, p_value, term_size, intersection_size, intersection.

For clusterProfiler GSEA outputs, this script uses core_enrichment as the
leading-edge signal gene set. The resulting interpretation is therefore
"leading-edge hitchhiking/redundancy" rather than ORA DEG-overlap redundancy.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


GENESET_OUTPUT_NAMES = {
    "hallmark",
    "CP_BIOCARTA",
    "CP_KEGG_MEDICUS",
    "CP_KEGG_LEGACY",
    "CP_REACTOME",
    "CP_WIKIPATHWAYS",
    "TFT_TFT_LEGACY",
    "TFT_GTRD",
    "GO_BP",
    "GO_CC",
    "GO_MF",
    "HPO",
    "C6",
    "IMMUNESIGDB",
}


@dataclass(frozen=True)
class GseaTask:
    data_type: str
    dataset_id: str
    plot_category: str
    analysis_name: str
    geneset_name: str
    gsea_file: Path
    output_dir: Path
    input_tsv: Path
    official_tsv: Path
    normalized_csv: Path
    stdout_file: Path
    stderr_file: Path
    gmt_file: Path


def sanitize_file_name(value: str, default: str = "analysis") -> str:
    value = str(value).strip()
    if not value:
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


def classify_gsea_path(gsea_file: Path, results_root: Path) -> tuple[str, str, str, str, str] | None:
    relative = gsea_file.relative_to(results_root)
    parts = relative.parts
    if len(parts) < 7 or parts[2] != "tables" or gsea_file.name != "gsea_result.csv":
        return None

    data_type = parts[0]
    dataset_id = parts[1]
    table_parts = parts[3:]

    try:
        gsea_index = table_parts.index("GSEA")
    except ValueError:
        return None

    if gsea_index + 2 >= len(table_parts):
        return None

    geneset_name = table_parts[gsea_index + 1]
    if geneset_name not in GENESET_OUTPUT_NAMES:
        return None

    prefix = table_parts[:gsea_index]
    if not prefix:
        return None

    if len(prefix) == 1:
        plot_category = "Main_DE"
        analysis_name = prefix[0]
    elif len(prefix) == 2 and prefix[0] == "ATF3_function":
        plot_category = "ATF3_correlation"
        analysis_name = prefix[1]
    elif len(prefix) == 3 and prefix[0] == "ATF3_function" and prefix[2] == "ATF3_high_low_DE":
        plot_category = "ATF3_high_low_DE"
        analysis_name = prefix[1]
    else:
        plot_category = sanitize_file_name("_".join(prefix[:-1]), default="Other_GSEA")
        analysis_name = prefix[-1]

    return data_type, dataset_id, plot_category, analysis_name, geneset_name


def discover_tasks(
    results_root: Path,
    gmt_dir: Path,
    datasets: set[str] | None,
    genesets: set[str] | None,
    categories: set[str] | None,
) -> list[GseaTask]:
    tasks: list[GseaTask] = []
    for gsea_file in sorted(results_root.glob("*/*/tables/**/GSEA/*/gsea_result.csv")):
        metadata = classify_gsea_path(gsea_file, results_root)
        if metadata is None:
            continue

        data_type, dataset_id, plot_category, analysis_name, geneset_name = metadata
        if datasets is not None and dataset_id not in datasets:
            continue
        if genesets is not None and geneset_name not in genesets:
            continue
        if categories is not None and plot_category not in categories:
            continue

        gmt_file = gmt_dir / f"{sanitize_file_name(geneset_name)}.gmt"
        output_dir = (
            results_root
            / data_type
            / dataset_id
            / "tables"
            / "PathwayDenester"
            / sanitize_file_name(plot_category)
            / sanitize_file_name(analysis_name)
            / sanitize_file_name(geneset_name)
        )

        tasks.append(
            GseaTask(
                data_type=data_type,
                dataset_id=dataset_id,
                plot_category=plot_category,
                analysis_name=analysis_name,
                geneset_name=geneset_name,
                gsea_file=gsea_file,
                output_dir=output_dir,
                input_tsv=output_dir / "pathway_denester_input.tsv",
                official_tsv=output_dir / "pathway_denester_official.tsv",
                normalized_csv=output_dir / "pathway_denester_result.csv",
                stdout_file=output_dir / "pathway_denester_stdout.log",
                stderr_file=output_dir / "pathway_denester_stderr.log",
                gmt_file=gmt_file,
            )
        )

    tasks

    return tasks


def ensure_official_tool(tool_script: Path, clone_dir: Path, clone_if_missing: bool) -> Path:
    if tool_script.exists():
        return tool_script

    candidate = clone_dir / "PathwayDenester.py"
    if candidate.exists():
        return candidate

    if not clone_if_missing:
        raise FileNotFoundError(
            f"PathwayDenester.py was not found at {tool_script} or {candidate}"
        )

    clone_dir.parent.mkdir(parents=True, exist_ok=True)
    if not clone_dir.exists():
        subprocess.run(
            ["git", "clone", "--depth", "1", "https://github.com/Helmy-Lab/PathwayDenester", str(clone_dir)],
            check=True,
        )

    if not candidate.exists():
        raise FileNotFoundError(f"Could not find cloned official tool: {candidate}")

    return candidate


def split_core_enrichment(value: object) -> list[str]:
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


def prepare_pathway_denester_input(task: GseaTask) -> tuple[int, int]:
    dat = pd.read_csv(task.gsea_file)
    required = {"ID", "Description", "pvalue", "core_enrichment"}
    missing = sorted(required - set(dat.columns))
    if missing:
        raise ValueError(f"Missing columns in {task.gsea_file}: {', '.join(missing)}")

    if "setSize" not in dat.columns:
        dat["setSize"] = pd.NA

    rows: list[dict[str, object]] = []
    for _, row in dat.iterrows():
        genes = split_core_enrichment(row.get("core_enrichment"))
        term_id = str(row.get("ID", "")).strip()
        term_name = str(row.get("Description", "")).strip()
        p_value = pd.to_numeric(row.get("pvalue"), errors="coerce")
        term_size = pd.to_numeric(row.get("setSize"), errors="coerce")

        if not term_id or not term_name or pd.isna(p_value) or not genes:
            continue

        rows.append(
            {
                "term_id": term_id,
                "term_name": term_name,
                "p_value": float(p_value),
                "term_size": "" if pd.isna(term_size) else int(term_size),
                "intersection_size": len(genes),
                "intersection": ",".join(genes),
            }
        )

    task.output_dir.mkdir(parents=True, exist_ok=True)
    prepared = pd.DataFrame(rows)
    prepared.to_csv(task.input_tsv, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    return len(dat), len(prepared)


def run_official_pathway_denester(
    task: GseaTask,
    official_script: Path,
    pval_threshold: float,
    to_test_threshold: float,
    top_n: int,
) -> None:
    cmd = [
        sys.executable,
        str(official_script),
        str(task.input_tsv),
        str(task.gmt_file),
        "--output_address",
        str(task.official_tsv),
        "--pval_threshold",
        str(pval_threshold),
        "--to_test_threshold",
        str(to_test_threshold),
        "--top_n",
        str(top_n),
    ]

    with task.stdout_file.open("w", encoding="utf-8") as stdout, task.stderr_file.open("w", encoding="utf-8") as stderr:
        subprocess.run(cmd, stdout=stdout, stderr=stderr, check=True)


def normalize_official_output(task: GseaTask) -> tuple[int, int, int, int]:
    dat = pd.read_csv(task.official_tsv, sep="\t")
    dat.insert(0, "Data_Type", task.data_type)
    dat.insert(1, "Dataset_ID", task.dataset_id)
    dat.insert(2, "Plot_Category", task.plot_category)
    dat.insert(3, "Analysis_Name", task.analysis_name)
    dat.insert(4, "GeneSet_Name", task.geneset_name)
    dat.insert(5, "GSEA_Result_File", str(task.gsea_file))

    dat.to_csv(task.normalized_csv, index=False)

    filtered = dat["Filtered"].astype(str).str.lower() if "Filtered" in dat.columns else pd.Series([], dtype=str)
    kept = int((filtered == "keep").sum())
    excluded = int((filtered == "exclude").sum())
    reciprocal = int(dat.get("Is Reciprocal", pd.Series([], dtype=str)).astype(str).str.lower().eq("true").sum())
    return len(dat), kept, excluded, reciprocal


def run_one_task(
    task: GseaTask,
    official_script: Path,
    pval_threshold: float,
    to_test_threshold: float,
    top_n: int,
    refresh: bool,
) -> dict[str, object]:
    base_record = {
        "Data_Type": task.data_type,
        "Dataset_ID": task.dataset_id,
        "Plot_Category": task.plot_category,
        "Analysis_Name": task.analysis_name,
        "GeneSet_Name": task.geneset_name,
        "GSEA_Result_File": str(task.gsea_file),
        "GMT_File": str(task.gmt_file),
        "PathwayDenester_Input_File": str(task.input_tsv),
        "PathwayDenester_Official_File": str(task.official_tsv),
        "PathwayDenester_Result_File": str(task.normalized_csv),
        "Stdout_File": str(task.stdout_file),
        "Stderr_File": str(task.stderr_file),
        "PValue_Threshold": pval_threshold,
        "To_Test_Threshold": to_test_threshold,
        "Top_N": top_n,
    }

    try:
        if not task.gmt_file.exists():
            raise FileNotFoundError(f"Missing GMT file: {task.gmt_file}")

        if task.normalized_csv.exists() and not refresh:
            dat = pd.read_csv(task.normalized_csv)
            filtered = dat["Filtered"].astype(str).str.lower() if "Filtered" in dat.columns else pd.Series([], dtype=str)
            evaluated_terms = len(dat)
            kept_terms = int((filtered == "keep").sum())
            excluded_terms = int((filtered == "exclude").sum())
            reciprocal_terms = int(dat.get("Is Reciprocal", pd.Series([], dtype=str)).astype(str).str.lower().eq("true").sum())
            input_terms = evaluated_terms
            prepared_terms = evaluated_terms
            status = "Skipped: existing result"
        else:
            input_terms, prepared_terms = prepare_pathway_denester_input(task)
            if prepared_terms < 2:
                raise ValueError("Fewer than two valid GSEA terms with core_enrichment.")
            run_official_pathway_denester(
                task=task,
                official_script=official_script,
                pval_threshold=pval_threshold,
                to_test_threshold=to_test_threshold,
                top_n=top_n,
            )
            evaluated_terms, kept_terms, excluded_terms, reciprocal_terms = normalize_official_output(task)
            status = "OK"

        excluded_percent = 100 * excluded_terms / evaluated_terms if evaluated_terms else 0.0
        return {
            **base_record,
            "Input_GSEA_Terms": input_terms,
            "Prepared_Terms": prepared_terms,
            "Evaluated_Terms": evaluated_terms,
            "Kept_Terms": kept_terms,
            "Excluded_Terms": excluded_terms,
            "Reciprocal_Terms": reciprocal_terms,
            "Excluded_Percent": round(excluded_percent, 4),
            "Status": status,
            "Error_Message": "",
        }
    except Exception as exc:  # keep batch jobs moving and summarize failures
        return {
            **base_record,
            "Input_GSEA_Terms": 0,
            "Prepared_Terms": 0,
            "Evaluated_Terms": 0,
            "Kept_Terms": 0,
            "Excluded_Terms": 0,
            "Reciprocal_Terms": 0,
            "Excluded_Percent": 0.0,
            "Status": "ERROR",
            "Error_Message": str(exc),
        }


def write_summary_tables(summary: pd.DataFrame, results_root: Path) -> None:
    global_dir = results_root / "PathwayDenester_summary"
    global_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(global_dir / "summary.csv", index=False)

    for (data_type, dataset_id), dat in summary.groupby(["Data_Type", "Dataset_ID"], sort=True):
        out_dir = results_root / data_type / dataset_id / "tables" / "PathwayDenester_summary"
        out_dir.mkdir(parents=True, exist_ok=True)
        dat.to_csv(out_dir / "summary.csv", index=False)


def run_tasks(args: argparse.Namespace) -> pd.DataFrame:
    results_root = Path(args.results_root)
    gmt_dir = Path(args.gmt_dir)
    official_script = ensure_official_tool(
        tool_script=Path(args.tool_script),
        clone_dir=Path(args.clone_dir),
        clone_if_missing=args.clone_if_missing,
    )

    tasks = discover_tasks(
        results_root=results_root,
        gmt_dir=gmt_dir,
        datasets=parse_csv_arg(args.datasets),
        genesets=parse_csv_arg(args.genesets),
        categories=parse_csv_arg(args.categories),
    )

    if not tasks:
        raise SystemExit("No GSEA result tables were found for PathwayDenester.")

    workers = max(int(args.workers), 1)
    print(f"PathwayDenester official script: {official_script}")
    print(f"Selected GSEA tables: {len(tasks)}")
    print(f"Workers: {workers}")
    print(f"pval_threshold: {args.pval_threshold}")
    print(f"to_test_threshold: {args.to_test_threshold}")

    records: list[dict[str, object]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_task = {
            executor.submit(
                run_one_task,
                task,
                official_script,
                args.pval_threshold,
                args.to_test_threshold,
                args.top_n,
                args.refresh,
            ): task
            for task in tasks
        }

        completed = 0
        for future in concurrent.futures.as_completed(future_to_task):
            completed += 1
            record = future.result()
            records.append(record)
            if completed == 1 or completed % max(1, min(20, len(tasks) // 10)) == 0 or completed == len(tasks):
                print(f"Completed {completed}/{len(tasks)}")

    summary = pd.DataFrame(records)
    summary = summary.sort_values(
        by=["Data_Type", "Dataset_ID", "Plot_Category", "Analysis_Name", "GeneSet_Name"]
    )
    write_summary_tables(summary, results_root)
    return summary


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch-run official PathwayDenester on all project GSEA results."
    )
    parser.add_argument("--results-root", default="results")
    parser.add_argument("--gmt-dir", default="data/reference/pathway_denester/msigdb_symbol_gmt")
    parser.add_argument("--tool-script", default="temporary/tools/PathwayDenester/PathwayDenester.py")
    parser.add_argument("--clone-dir", default="temporary/tools/PathwayDenester")
    parser.add_argument("--clone-if-missing", action="store_true")
    parser.add_argument("--datasets", default="all", help="Comma-separated dataset IDs or all.")
    parser.add_argument("--genesets", default="all", help="Comma-separated geneset output names or all.")
    parser.add_argument("--categories", default="all", help="Comma-separated plot categories or all.")
    parser.add_argument("--workers", type=int, default=max(1, min((os.cpu_count() or 2) // 2, 6)))
    parser.add_argument("--pval-threshold", type=float, default=0.05)
    parser.add_argument("--to-test-threshold", type=float, default=0.5)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--refresh", action="store_true")
    return parser.parse_args(list(argv))


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    summary = run_tasks(args)
    ok_count = int(summary["Status"].astype(str).eq("OK").sum())
    skipped_count = int(summary["Status"].astype(str).str.startswith("Skipped").sum())
    error_count = int(summary["Status"].astype(str).eq("ERROR").sum())
    print("\nPathwayDenester batch finished.")
    print(f"OK: {ok_count}")
    print(f"Skipped: {skipped_count}")
    print(f"Errors: {error_count}")
    print(f"Global summary: {Path(args.results_root) / 'PathwayDenester_summary' / 'summary.csv'}")
    return 1 if error_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
