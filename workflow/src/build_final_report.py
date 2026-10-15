#!/usr/bin/env python3
"""
Build a scrollable HTML report and a candidate-covariation TSV from a
DarwinsRNAHunt motif-discovery run.

For each cluster this pulls together:
  - the R-scape .cov file (ranked covarying base pairs)
  - the captured R-scape stdout (for the power / summary text)
  - the R2R CaCoFold structure diagram (PDF), rendered inline as a PNG
  - cmstat output for the combined calibrated covariance models

Only standard library is used, plus an optional shell-out to `cmstat`
(Infernal) and `pdftoppm` (poppler-utils) if they're on PATH. Both are
soft dependencies: if either is missing, the report is still built, just
without that section (cmstat text becomes an explanatory note; R2R
diagrams become a download link instead of an inline image).

Usage:
    python build_final_report.py \\
        --rscape-dirs results/.../03_structure/cluster_5_rscape_output \\
                      results/.../03_structure/cluster_12_rscape_output \\
        --cm-stats results/.../04_cms/5_calibrated.cm.stats \\
                   results/.../04_cms/12_calibrated.cm.stats \\
        --target-domain 16S --region upstream \\
        --html-out results/.../final_report.html \\
        --tsv-out results/.../candidate_covariations.tsv

Each cluster's ID is parsed from its rscape-output directory name; the
matching struc.txt is found as a sibling of that directory (both live
under .../03_structure/). --cm-stats files are matched to clusters by
their own filename (<id>_calibrated.cm.stats), not by list position, so
the two lists don't need to be given in the same order. Nothing is
discovered by globbing a base directory - the caller decides exactly
which clusters go into the report.
"""

import argparse
import base64
import csv
import html
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path


# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------

@dataclass
class CovPair:
    in_structure: str
    left: int
    right: int
    score: float
    evalue: float
    subs: float
    power: float


@dataclass
class ClusterResult:
    cluster_id: str
    cov_pairs: list = field(default_factory=list)
    cov_summary: dict = field(default_factory=dict)
    cov_raw_comments: list = field(default_factory=list)
    struc_text: str = ""
    r2r_pdf: Path = None
    r2r_png_b64: str = None
    cmstat_text: str = ""

    @property
    def n_pairs(self):
        return len(self.cov_pairs)

    @property
    def best_evalue(self):
        if not self.cov_pairs:
            return None
        return min(p.evalue for p in self.cov_pairs)


# --------------------------------------------------------------------------
# Parsing R-scape .cov files
# --------------------------------------------------------------------------
# Column layout per the R-scape user guide: flag, left, right, score,
# E-value, [substitutions, power]. Two "#" comment lines per alignment
# carry nseq/alen/avgid/nbpairs and the covariation-method/threshold
# summary - format varies a bit between R-scape versions, so we pull
# what we can with regexes and always keep the raw text too.

_SUMMARY_PATTERNS = {
    "nseq": r"nseq\s+(\d+)",
    "alen": r"alen\s+(\d+)",
    "avgid": r"avgid\s+([\d.]+)",
    "nbpairs": r"nbpairs\s+(\d+)",
    "evalue_threshold": r"E-?value\s*(?:target|threshold)?\D{0,10}?([\d.]+)",
    "n_cov_nonbp": r"covarying non[- ]?base\s*pairs\D*(\d+)",
    "n_cov_bp": r"covarying base\s*pairs\D*(\d+)",
}


def _extract_summary_fields(comment_line: str, summary: dict) -> None:
    for key, pattern in _SUMMARY_PATTERNS.items():
        if key in summary:
            continue
        match = re.search(pattern, comment_line, re.IGNORECASE)
        if match:
            summary[key] = match.group(1)


def parse_cov_file(path: Path):
    """ example .cov file:
    #
    # Method Target_E-val [cov_min,cov_max] [FP | TP True Found | Sen PPV F] 
    # GTp    0.05           [-9.95,160.57]    [0 | 2 98 2 | 2.04 100.00 4.00] 
    #
    # in_CaCoFold in_given  left_pos      right_pos  score            E-value       pvalue       substitutions    power
    #------------------------------------------------------------------------------------------------------------------------
    *	~	      7058	      7098	130.07956	0.00923598	2.17389e-07	29		0.73
    *	~	      7105	      7121	160.56991	0.000276359	6.50471e-09	24		0.66
    """

    pairs, summary, raw_comments = [], {}, []
    if not path or not path.exists() or path.stat().st_size == 0:
        return pairs, summary, raw_comments

    for line in path.read_text().splitlines():
        line = line.rstrip()
        if not line:
            continue
        if line.startswith("#"):
            raw_comments.append(line)
            _extract_summary_fields(line, summary)
            continue
        fields = line.split()
        if len(fields) < 5:
            continue
        try:
            in_cacofold = fields[0]
            in_given = fields[1]
            left = int(fields[2])
            right = int(fields[3])
            score = float(fields[4])
            evalue = float(fields[5])
            pvalue = float(fields[6])
            subs = float(fields[7])
            power = float(fields[8])
        except ValueError:
            continue
        pairs.append(CovPair(in_cacofold, left, right, score, evalue, subs, power))

    return pairs, summary, raw_comments


# --------------------------------------------------------------------------
# Locating per-cluster files
# --------------------------------------------------------------------------
# The caller (Snakemake) passes one R-scape output directory per cluster
# explicitly - no filesystem discovery/globbing for *which* clusters exist.
# Everything else (the .cov file, the R2R pdf, the struc.txt) is located
# relative to that directory, since they all live alongside it as siblings
# under 03_structure/.

_RSCAPE_DIR_RE = re.compile(r"cluster_(.+)_rscape_output$")


def cluster_id_from_rscape_dir(rscape_dir: Path) -> str:
    m = _RSCAPE_DIR_RE.match(rscape_dir.name)
    if not m:
        raise ValueError(
            f"Can't extract a cluster ID from directory name {rscape_dir.name!r}; "
            "expected it to look like 'cluster_<id>_rscape_output'."
        )
    return m.group(1)


def find_struc_text(rscape_dir: Path, cluster_id: str) -> str:
    # struc.txt is a sibling of the rscape_output dir, both directly under 03_structure/
    path = rscape_dir.parent / f"cluster_{cluster_id}_struc.txt"
    return path.read_text() if path.exists() else ""


def find_cov_file(rscape_dir: Path, cluster_id: str):
    candidates = sorted(rscape_dir.glob(f"cluster_{cluster_id}_1.cacofold.cov"))
    return candidates[0] if candidates else None


def find_r2r_pdf(rscape_dir: Path, cluster_id: str):
    candidates = sorted(rscape_dir.glob(f"cluster_{cluster_id}_*.cacofold.R2R.sto.pdf"))
    return candidates[0] if candidates else None


def render_pdf_to_png_b64(pdf_path: Path, dpi: int = 150):
    """Render page 1 of the R2R PDF to a base64 PNG so the HTML report is a
    single self-contained file. Returns None (caller falls back to a link)
    if pdftoppm isn't available or conversion fails."""
    if pdf_path is None or shutil.which("pdftoppm") is None:
        return None
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = Path(tmp) / "page"
        try:
            subprocess.run(
                ["pdftoppm", "-png", "-r", str(dpi), "-singlefile", str(pdf_path), str(out_prefix)],
                check=True,
            )
        except subprocess.CalledProcessError:
            return None
        png_path = out_prefix.with_suffix(".png")
        if not png_path.exists():
            return None
        return base64.b64encode(png_path.read_bytes()).decode("ascii")


_CM_STATS_RE = re.compile(r"^(.+)_calibrated\.cm\.stats$")


def cluster_id_from_cm_stats_path(path: Path) -> str:
    m = _CM_STATS_RE.match(path.name)
    if not m:
        raise ValueError(
            f"Can't extract a cluster ID from cm-stats filename {path.name!r}; "
            "expected it to look like '<id>_calibrated.cm.stats'."
        )
    return m.group(1)


def read_cm_stats(cm_stats_paths) -> dict:
    """Read the cmstat output already captured by calibrate_cm_model, keyed
    by cluster ID. No cmstat subprocess needed here - it was already run
    once upstream, so we just display what's there."""
    by_cluster = {}
    for path in cm_stats_paths:
        cid = cluster_id_from_cm_stats_path(path)
        if not path.exists() or path.stat().st_size == 0:
            by_cluster[cid] = "(no cmstat output - calibration may have been skipped for this cluster)"
            continue
        by_cluster[cid] = path.read_text()
    return by_cluster


def gather_cluster_results(rscape_dirs, cm_stats_by_cluster, dpi: int):
    results = []
    for rscape_dir in rscape_dirs:
        cid = cluster_id_from_rscape_dir(rscape_dir)
        if not rscape_dir.exists():
            print(f"warning: {rscape_dir} does not exist, skipping cluster {cid}", file=sys.stderr)
            continue

        cov_path = find_cov_file(rscape_dir, cid)
        pairs, summary, raw_comments = parse_cov_file(cov_path)

        r2r_pdf = find_r2r_pdf(rscape_dir, cid)
        r2r_png_b64 = render_pdf_to_png_b64(r2r_pdf, dpi)

        results.append(ClusterResult(
            cluster_id=cid,
            cov_pairs=pairs,
            cov_summary=summary,
            cov_raw_comments=raw_comments,
            struc_text=find_struc_text(rscape_dir, cid),
            r2r_pdf=r2r_pdf,
            r2r_png_b64=r2r_png_b64,
            cmstat_text=cm_stats_by_cluster.get(cid, "(no cmstat file provided for this cluster)"),
        ))

    # Most promising first: most covarying pairs, then strongest E-value.
    results.sort(key=lambda r: (-r.n_pairs, r.best_evalue if r.best_evalue is not None else float("inf")))
    return results


# --------------------------------------------------------------------------
# TSV output - only clusters with at least one significant covarying pair
# --------------------------------------------------------------------------

def write_tsv(results, target_domain, region, tsv_out: Path):
    tsv_out.parent.mkdir(parents=True, exist_ok=True)
    with tsv_out.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "target_domain", "region", "cluster_id", "in_structure",
            "left_pos", "right_pos", "score", "evalue", "substitutions", "power",
        ])
        for r in results:
            for p in r.cov_pairs:
                writer.writerow([
                    target_domain, region, r.cluster_id, p.in_structure,
                    p.left, p.right, p.score, p.evalue, p.subs, p.power,
                ])


# --------------------------------------------------------------------------
# HTML output
# --------------------------------------------------------------------------

_CSS = """
body { font-family: -apple-system, Helvetica, Arial, sans-serif; margin: 2em auto; max-width: 960px; color: #222; }
h1 { border-bottom: 2px solid #333; padding-bottom: 0.2em; }
h2 { margin-top: 2.5em; border-top: 1px solid #ccc; padding-top: 1em; }
table { border-collapse: collapse; width: 100%; margin: 0.5em 0 1em; font-size: 0.92em; }
th, td { border: 1px solid #ddd; padding: 4px 8px; text-align: right; }
th { background: #f2f2f2; text-align: center; }
td:first-child, th:first-child { text-align: left; }
tr.promising { background: #eaf7ea; }
.pill { display: inline-block; padding: 2px 8px; border-radius: 10px; font-size: 0.85em; margin-left: 0.5em; }
.pill-good { background: #cdeecd; }
.pill-none { background: #eee; color: #777; }
img.r2r { max-width: 100%; border: 1px solid #ccc; }
pre { background: #f7f7f7; padding: 0.8em; overflow-x: auto; font-size: 0.85em; }
a.jump { font-size: 0.9em; }
.overview td, .overview th { text-align: center; }
.overview td:first-child { text-align: left; }
"""


def _fmt(value, digits=3):
    if value is None:
        return "-"
    try:
        f = float(value)
    except (TypeError, ValueError):
        return html.escape(str(value))
    if f != f:  # NaN
        return "-"
    return f"{f:.{digits}g}"


def build_overview_table(results):
    rows = ['<table class="overview"><tr><th>Cluster</th><th># covarying pairs</th>'
            '<th>Best E-value</th><th>nseq</th><th>alen</th><th>avgid</th><th></th></tr>']
    for r in results:
        promising = r.n_pairs > 0
        cls = ' class="promising"' if promising else ""
        pill = '<span class="pill pill-good">candidate</span>' if promising else '<span class="pill pill-none">no covariation</span>'
        rows.append(
            f'<tr{cls}><td><a class="jump" href="#cluster-{html.escape(r.cluster_id)}">'
            f'{html.escape(r.cluster_id)}</a></td>'
            f'<td>{r.n_pairs}</td><td>{_fmt(r.best_evalue)}</td>'
            f'<td>{html.escape(r.cov_summary.get("nseq", "-"))}</td>'
            f'<td>{html.escape(r.cov_summary.get("alen", "-"))}</td>'
            f'<td>{html.escape(r.cov_summary.get("avgid", "-"))}</td>'
            f'<td>{pill}</td></tr>'
        )
    rows.append("</table>")
    return "\n".join(rows)


def build_cluster_section(r: ClusterResult):
    parts = [f'<h2 id="cluster-{html.escape(r.cluster_id)}">Cluster {html.escape(r.cluster_id)}</h2>']

    summary = r.cov_summary
    if summary:
        items = "".join(f"<li><b>{html.escape(k)}</b>: {html.escape(str(v))}</li>" for k, v in summary.items())
        parts.append(f"<ul>{items}</ul>")
    else:
        parts.append("<p><i>No .cov summary line could be parsed for this cluster.</i></p>")

    if r.cov_pairs:
        rows_sorted = sorted(r.cov_pairs, key=lambda p: p.evalue)
        table = ['<table><tr><th>In structure</th><th>i</th><th>j</th><th>Score</th>'
                 '<th>E-value</th><th>Substitutions</th><th>Power</th></tr>']
        for p in rows_sorted:
            table.append(
                f"<tr><td>{html.escape(p.in_structure)}</td><td>{p.left}</td><td>{p.right}</td>"
                f"<td>{_fmt(p.score)}</td><td>{_fmt(p.evalue, 3)}</td>"
                f"<td>{_fmt(p.subs)}</td><td>{_fmt(p.power)}</td></tr>"
            )
        table.append("</table>")
        parts.append("\n".join(table))
    else:
        parts.append("<p><i>No significant covarying pairs reported for this cluster.</i></p>")

    parts.append(f"<p><b>cmstat</b></p><pre>{html.escape(r.cmstat_text)}</pre>")

    if r.r2r_png_b64:
        parts.append(
            f'<p><img class="r2r" alt="R2R structure diagram for cluster {html.escape(r.cluster_id)}" '
            f'src="data:image/png;base64,{r.r2r_png_b64}"></p>'
        )
    elif r.r2r_pdf:
        parts.append(f'<p><a href="{html.escape(str(r.r2r_pdf))}">R2R structure diagram (PDF)</a> '
                      f'- install poppler-utils (pdftoppm) to embed this inline.</p>')
    else:
        parts.append("<p><i>No R2R structure diagram found for this cluster.</i></p>")

    if r.struc_text.strip():
        parts.append("<details><summary>R-scape run output</summary>"
                      f"<pre>{html.escape(r.struc_text)}</pre></details>")

    if r.cov_raw_comments:
        parts.append("<details><summary>Raw .cov header lines</summary>"
                      f"<pre>{html.escape(chr(10).join(r.cov_raw_comments))}</pre></details>")

    return "\n".join(parts)


def write_html(results, target_domain, region, html_out: Path):
    html_out.parent.mkdir(parents=True, exist_ok=True)
    n_candidates = sum(1 for r in results if r.n_pairs > 0)

    body = [
        f"<h1>Motif discovery report: {html.escape(target_domain)} / {html.escape(region)}</h1>",
        f"<p>{len(results)} clusters analysed, {n_candidates} with at least one significant "
        "covarying pair (highlighted below). Clusters are ranked by number of covarying pairs, "
        "then by best E-value.</p>",
        build_overview_table(results),
    ]
    for r in results:
        body.append(build_cluster_section(r))

    doc = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Motif discovery report: {html.escape(target_domain)} / {html.escape(region)}</title>
<style>{_CSS}</style>
</head>
<body>
{chr(10).join(body)}
</body>
</html>
"""
    html_out.write_text(doc)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rscape-dirs", required=True, nargs="*", type=Path,
                         help="One or more 'cluster_<id>_rscape_output' directories, one per cluster "
                              "(e.g. results/.../03_structure/cluster_5_rscape_output). Can be empty "
                              "if no clusters passed upstream filtering - an empty report is written.")
    parser.add_argument("--cm-stats", nargs="*", type=Path, default=[],
                         help="Per-cluster '<id>_calibrated.cm.stats' files, as already produced by the "
                              "calibrate_cm_model rule. Matched to clusters by filename, not position.")
    parser.add_argument("--target-domain", default="")
    parser.add_argument("--region", default="")
    parser.add_argument("--html-out", required=True, type=Path)
    parser.add_argument("--tsv-out", required=True, type=Path)
    parser.add_argument("--dpi", type=int, default=150, help="Resolution for R2R PDF -> PNG rendering")
    args = parser.parse_args()

    if not args.rscape_dirs:
        print("No --rscape-dirs given - writing an empty report "
              "(check upstream clustering if this is unexpected).", file=sys.stderr)

    cm_stats_by_cluster = read_cm_stats(args.cm_stats)
    results = gather_cluster_results(args.rscape_dirs, cm_stats_by_cluster, args.dpi)

    write_tsv(results, args.target_domain, args.region, args.tsv_out)
    write_html(results, args.target_domain, args.region, args.html_out)

    n_candidates = sum(1 for r in results if r.n_pairs > 0)
    print(f"Wrote {args.html_out} and {args.tsv_out}: "
          f"{n_candidates}/{len(results)} clusters have significant covariation.")


if __name__ == "__main__":
    main()