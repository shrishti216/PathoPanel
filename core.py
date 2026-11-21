"""PathoPanel analysis engine (no Streamlit code in here, so it can be tested alone).

Each oligo is aligned end to end against both strands of every sequence, with no gaps
inside the oligo, and the placement with the fewest mismatches is kept. Mismatches are
then described by position, pair type and distance from the 3' end.
"""
from __future__ import annotations

import base64
import datetime
import hashlib
import io
import random
import re
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from Bio import SeqIO
from Bio.Align import PairwiseAligner

LABEL = {"Fwd": "Forward primer", "Rev": "Reverse primer", "Probe": "Probe"}
IS_PRIMER = {"Fwd": True, "Rev": True, "Probe": False}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
# Pair types from Stadhouders et al. 2010, J Mol Diagn 12(1):109-117 (primer/template pairs).
# That study examined the last bases of the primer; applying the weights to internal
# mismatches is an extension. Pairs not ranked in that paper are treated as intermediate.
MILD, SEVERE = {"AC", "GT"}, {"AA", "AG", "CC"}
CLASS_ORDER = ["Exact", "Minor", "Moderate", "High risk", "No site"]
STATUS_ORDER = CLASS_ORDER + ["No amplicon"]
COLORS = {"Exact": "#0F766E", "Minor": "#7FBF9E", "Moderate": "#E3A72F",
          "High risk": "#D9534F", "No site": "#9AA3A8", "No amplicon": "#5B6670"}

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -10
aligner.extend_gap_score = -10
aligner.end_deletion_score = 0      # the long sequence may overhang freely


# ------------------------------------------------------------------ input
def clean_oligo(text: str) -> str:
    return re.sub(r"\s+", "", text or "").upper().replace("U", "T")


def valid_oligo(seq: str) -> bool:
    return bool(seq) and re.fullmatch(r"[ACGT]+", seq) is not None


def serotype(text: str) -> str:
    m = re.search(r"(?:dengue virus|DENV)[\s_/-]*(?:type)?[\s_/-]*([1-4])", text, re.I)
    return "DENV" + m.group(1) if m else "unknown"


def read_fasta(raw: bytes):
    return list(SeqIO.parse(io.StringIO(raw.decode("utf-8", errors="ignore")), "fasta"))


def filter_and_sample(records, min_length=0, sample_size=0, balance=False, seed=42):
    n_in = len(records)
    if min_length:
        records = [r for r in records if len(r.seq) >= min_length]
    n_short = n_in - len(records)
    if sample_size and sample_size < len(records):
        rng = random.Random(seed)
        if balance:
            groups = {}
            for r in records:
                groups.setdefault(serotype(r.description), []).append(r)
            per = max(1, sample_size // len(groups))
            records = [r for g in groups.values() for r in rng.sample(g, min(per, len(g)))]
        else:
            records = rng.sample(records, sample_size)
    return records, n_in, n_short


# ----------------------------------------------------------------- scoring
def mismatch_weight(oligo_base, site_base, dist_from_3p, is_primer):
    template_base = COMP.get(site_base, "N")
    pair = "".join(sorted(oligo_base + template_base))
    type_w = 0.5 if pair in MILD else 2.0 if pair in SEVERE else 1.0
    if not is_primer:
        return type_w
    pos_w = 4.0 if dist_from_3p == 1 else 3.0 if dist_from_3p <= 3 else 2.0 if dist_from_3p <= 5 else 1.0
    return type_w * pos_w


def classify(n_mm, score):
    """Heuristic labels. The four-mismatch cutoff is an arbitrary threshold."""
    if n_mm == 0:
        return "Exact"
    if n_mm >= 4 or score > 5:
        return "High risk"
    return "Minor" if score <= 1.5 else "Moderate"


NO_SITE = {"class": "No site", "mm": [], "n_mm": None, "score": None, "term3": None,
           "strand": "", "start": None, "end": None, "site": "", "types": ""}


def find_site(key, oligo, seq, max_mm, max_unreadable):
    L = len(oligo)
    best = None
    for strand, target in (("+", seq), ("-", seq.reverse_complement())):
        try:
            s = int(aligner.align(target, oligo)[0].aligned[0][0][0])
        except (IndexError, ValueError):
            continue
        site = str(target[s:s + L])
        if len(site) < L:
            continue
        mm = [i for i in range(L) if site[i] in "ACGT" and site[i] != oligo[i]]
        if best is None or len(mm) < len(best["mm"]):
            start = s + 1 if strand == "+" else len(seq) - (s + L) + 1
            best = {"strand": strand, "site": site, "mm": mm, "start": start}
    if best is None or len(best["mm"]) > max_mm:
        return dict(NO_SITE)
    site, mm = best["site"], best["mm"]
    if sum(1 for ch in site if ch not in "ACGT") > max_unreadable:
        return dict(NO_SITE)            # too many unreadable bases to call the site
    score = sum(mismatch_weight(oligo[i], site[i], L - i, IS_PRIMER[key]) for i in mm)
    best.update({
        "n_mm": len(mm), "score": round(score, 2),
        "term3": sum(1 for i in mm if L - i <= 5),
        "class": classify(len(mm), score), "end": best["start"] + L - 1,
        "types": "; ".join(f"{i + 1}:{oligo[i]}/{COMP.get(site[i], 'N')}" for i in mm)})
    return best


def amplicon(hits, max_amp):
    f, r = hits.get("Fwd"), hits.get("Rev")
    if not f or not r or f["class"] == "No site" or r["class"] == "No site":
        return None, None
    if f["strand"] == r["strand"]:
        return None, False
    plus, minus = (f, r) if f["strand"] == "+" else (r, f)
    length = minus["end"] - plus["start"] + 1
    return length, 0 < length <= max_amp


def analyse(records, oligos, max_mm=6, max_amp=2000, max_unreadable=2, progress=None):
    rows = []
    site_log = {k: [] for k in oligos}
    total = len(records)
    for n, rec in enumerate(records, 1):
        seq = rec.seq.upper()
        hits = {k: find_site(k, o, seq, max_mm, max_unreadable) for k, o in oligos.items()}
        status = max((h["class"] for h in hits.values()), key=CLASS_ORDER.index)
        amp_len, amp_ok = amplicon(hits, max_amp)
        if amp_ok is False and status in ("Exact", "Minor"):
            status = "No amplicon"
        sero = serotype(rec.description)
        row = {"Sequence": rec.id, "Serotype": sero, "Length": len(seq), "Status": status,
               "Amplicon length": amp_len}
        for k, h in hits.items():
            row[f"{k} class"] = h["class"]
            row[f"{k} mismatches"] = h["n_mm"]
            row[f"{k} 3' end mismatches"] = h["term3"]
            row[f"{k} risk score"] = h["score"]
            row[f"{k} positions"] = ", ".join(str(i + 1) for i in h["mm"])
            row[f"{k} mismatch types"] = h["types"]
            row[f"{k} site start"] = h["start"]
            if h["class"] != "No site":
                site_log[k].append((sero, h["site"], h["mm"]))
        rows.append(row)
        if progress and (n % 50 == 0 or n == total):
            progress(n / total)
    return pd.DataFrame(rows), site_log


# ---------------------------------------------------------------- summaries
def pct(mask):
    return round(100 * mask.mean(), 1)


def serotype_order(df):
    known = [s for s in ["DENV1", "DENV2", "DENV3", "DENV4", "unknown"] if s in set(df["Serotype"])]
    return known + sorted(set(df["Serotype"]) - set(known))


def primer_keys(oligos):
    return [k for k in ("Fwd", "Rev") if k in oligos]


RULE_LABELS = [
    "A. Every primer has a binding site (any number of mismatches)",
    "B. A, and the last 5 bases at each 3' end match",
    "C. B, and at most 4 mismatches in each primer",
    "D. B, and at most 2 mismatches in each primer",
]


def rules_table(df, oligos):
    pk = primer_keys(oligos)
    if not pk:
        return None
    has_site = np.logical_and.reduce([df[f"{k} mismatches"].notna() for k in pk])
    end_ok = np.logical_and.reduce([df[f"{k} 3' end mismatches"] == 0 for k in pk])
    le = lambda n: np.logical_and.reduce([df[f"{k} mismatches"] <= n for k in pk])
    masks = [has_site, end_ok, end_ok & le(4), end_ok & le(2)]
    table = pd.DataFrame({lab: df.assign(ok=m).groupby("Serotype")["ok"].mean().mul(100).round(1)
                          for lab, m in zip(RULE_LABELS, masks)}).T
    table = table[[c for c in serotype_order(df) if c in table.columns]]
    table["All sequences"] = [round(100 * m.mean(), 1) for m in masks]
    return table


def scorecard(df, oligos):
    def one(d):
        out = {"Sequences": len(d)}
        for k in oligos:
            out[f"{LABEL[k]}: exact (%)"] = pct(d[f"{k} class"] == "Exact")
            out[f"{LABEL[k]}: exact or minor (%)"] = pct(d[f"{k} class"].isin(["Exact", "Minor"]))
        out["Whole assay: exact or minor (%)"] = pct(d["Status"].isin(["Exact", "Minor"]))
        return pd.Series(out)
    card = df.groupby("Serotype").apply(one, include_groups=False)
    card = card.loc[serotype_order(df)]
    card.loc["All sequences"] = one(df)
    card_t = card.T.astype(object)
    card_t.loc["Sequences"] = card_t.loc["Sequences"].astype(int)
    return card_t


def three_prime_summary(df, oligos):
    out = {}
    for k in primer_keys(oligos):
        col = f"{k} 3' end mismatches"
        has = df[df[col] > 0]
        by = pd.DataFrame({"Sequences with a 3' end mismatch": has["Serotype"].value_counts()})
        totals = df["Serotype"].value_counts()
        by["Share of serotype (%)"] = (100 * by.iloc[:, 0] / totals.reindex(by.index)).round(1)
        by = by.reindex([s for s in serotype_order(df) if s in by.index])
        out[k] = {"n": len(has), "by": by,
                  "positions": has[f"{k} positions"].value_counts().head(8).rename("Sequences").to_frame(),
                  "types": has[f"{k} mismatch types"].value_counts().head(8).rename("Sequences").to_frame()}
    return out


def variant_table(site_log, key, top=8):
    if not site_log[key]:
        return None
    vdf = pd.DataFrame(site_log[key], columns=["Serotype", "Binding site", "mm"])
    vdf["Mismatches"] = vdf["mm"].str.len()
    top_t = (vdf.groupby("Binding site").agg(Sequences=("Serotype", "size"), Mismatches=("Mismatches", "first"))
             .sort_values("Sequences", ascending=False).head(top))
    top_t["Share (%)"] = (100 * top_t["Sequences"] / len(vdf)).round(1)
    return top_t.join(pd.crosstab(vdf["Binding site"], vdf["Serotype"]))


IUPAC_FROM_SET = {frozenset(k): v for k, v in {
    "A": "A", "C": "C", "G": "G", "T": "T", "AG": "R", "CT": "Y", "CG": "S", "AT": "W",
    "GT": "K", "AC": "M", "CGT": "B", "AGT": "D", "ACT": "H", "ACG": "V", "ACGT": "N"}.items()}


def redesign(oligos, site_log, threshold):
    rows = []
    for k, oligo in oligos.items():
        if not site_log[k]:
            continue
        L = len(oligo)
        sites = Counter(s for _, s, _ in site_log[k])
        total = sum(sites.values())
        allowed = []
        for i in range(L):
            counts = Counter()
            for site, c in sites.items():
                if site[i] in "ACGT":
                    counts[site[i]] += c
            allowed.append({oligo[i]} | {b for b, c in counts.items() if c / total >= threshold})
        changed = [i + 1 for i in range(L) if allowed[i] != {oligo[i]}]
        now = sum(c for s, c in sites.items() if s == oligo)
        after = sum(c for s, c in sites.items() if all(s[i] in allowed[i] for i in range(L)))
        fold = int(np.prod([len(b) for b in allowed]))
        caution = []
        if fold > 16:
            caution.append("very high degeneracy")
        if any(L - (p - 1) <= 3 for p in changed):
            caution.append("change within the last 3 bases at the 3' end")
        rows.append({"Oligo": LABEL[k], "Current": oligo,
                     "Suggested": "".join(IUPAC_FROM_SET[frozenset(b)] for b in allowed),
                     "Changed positions": ", ".join(map(str, changed)) or "none",
                     "Degeneracy": f"{fold}-fold",
                     "Exact sites now (%)": round(100 * now / total, 1),
                     "Exact sites with suggestion (%)": round(100 * after / total, 1),
                     "Caution": "; ".join(caution)})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ figures
def fig_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    return base64.b64encode(buf.getvalue()).decode()


def fig_status(df, oligos):
    sero = serotype_order(df)
    panels = [(f"{k} class", LABEL[k]) for k in oligos] + [("Status", "Whole assay")]
    fig, axes = plt.subplots(len(panels), 1, figsize=(9, 2.1 * len(panels)), squeeze=False)
    for ax, (col, title) in zip(axes[:, 0], panels):
        left = np.zeros(len(sero))
        for status in STATUS_ORDER:
            vals = np.array([100 * (df.loc[df["Serotype"] == s, col] == status).mean() for s in sero])
            ax.barh(range(len(sero)), vals, left=left, color=COLORS[status],
                    label=status if col == "Status" else None)
            left += vals
        ax.set_yticks(range(len(sero)))
        ax.set_yticklabels([f"{s} (n={int((df['Serotype'] == s).sum()):,})" for s in sero])
        ax.invert_yaxis(); ax.set_xlim(0, 100); ax.set_title(title, loc="left", fontsize=11)
        ax.set_xlabel("Share of sequences (%)" if col == "Status" else "")
    axes[-1, 0].legend(ncol=6, loc="upper center", bbox_to_anchor=(0.5, -0.5), frameon=False, fontsize=8)
    plt.tight_layout()
    return fig


def fig_landscape(df, site_log, oligos):
    sero = serotype_order(df)
    fig, axes = plt.subplots(len(oligos), 1, figsize=(10, 2.0 + 1.6 * len(oligos)), squeeze=False)
    for ax, (k, oligo) in zip(axes[:, 0], oligos.items()):
        L = len(oligo)
        grid = np.zeros((len(sero), L))
        for r_i, s in enumerate(sero):
            found = [mm for sr, _, mm in site_log[k] if sr == s]
            for mm in found:
                for i in mm:
                    grid[r_i, i] += 1
            if found:
                grid[r_i] = 100 * grid[r_i] / len(found)
        im = ax.imshow(grid, cmap="YlOrRd", vmin=0, vmax=100, aspect="auto")
        ax.set_xticks(range(L)); ax.set_xticklabels(list(oligo), fontsize=7)
        ax.set_yticks(range(len(sero))); ax.set_yticklabels(sero)
        for r_i in range(len(sero)):
            for c_i in range(L):
                if grid[r_i, c_i] >= 5:
                    ax.text(c_i, r_i, f"{grid[r_i, c_i]:.0f}", ha="center", va="center", fontsize=6)
        if IS_PRIMER[k]:
            ax.add_patch(Rectangle((L - 5.5, -0.5), 5, len(sero), fill=False, ec="black", lw=1.6))
        ax.set_title(f"{LABEL[k]}: % of sequences differing at each base (5' on the left)",
                     loc="left", fontsize=10)
        plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    plt.tight_layout()
    return fig


def fig_amplicon(df, expected):
    if not df["Amplicon length"].notna().any():
        return None
    sero = serotype_order(df)
    fig, ax = plt.subplots(figsize=(9, 3))
    data = [df.loc[(df["Serotype"] == s) & df["Amplicon length"].notna(), "Amplicon length"] for s in sero]
    ax.hist(data, bins=40, stacked=True, label=sero)
    if expected:
        ax.axvline(expected, color="black", ls="--", lw=1, label=f"expected {expected} bp")
    ax.set_xlabel("Amplicon length (bp)"); ax.set_ylabel("Sequences"); ax.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    return fig


# ------------------------------------------------------------------- report
def build_report(meta, oligos, df, rules, card, repair, figs):
    img = lambda n, alt: (f'<img alt="{alt}" style="max-width:100%" src="data:image/png;base64,{figs[n]}">'
                          if figs.get(n) else "")
    oligo_rows = "".join(f"<tr><td>{LABEL[k]}</td><td><code>{v}</code></td><td>{len(v)} nt</td></tr>"
                         for k, v in oligos.items())
    allc = rules["All sequences"] if rules is not None else None
    summary = ""
    if allc is not None:
        summary = (f"<div class='note'><b>{allc.iloc[0]}%</b> of {len(df):,} sequences carry a binding site for "
                   f"every primer, and <b>{allc.iloc[1]}%</b> also match in the last five bases at each 3' end. "
                   f"Compatibility ranges from {allc.iloc[3]}% to {allc.iloc[0]}% depending on the mismatch "
                   f"tolerance assumed. This is sequence compatibility only. A match does not prove amplification, "
                   f"and a mismatch does not prove diagnostic failure.</div>")
    sample_note = f" (random sample of {meta['sample_size']}, seed 42)" if meta["sample_size"] else ""
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>PathoPanel report: {meta['assay']}</title>
<style>
body{{font-family:Segoe UI,Helvetica,Arial,sans-serif;color:#14232B;max-width:900px;margin:2rem auto;padding:0 1rem;line-height:1.5}}
h2{{margin-top:2rem;border-bottom:2px solid #0F766E;padding-bottom:.2rem}}
table{{border-collapse:collapse;margin:.6rem 0;font-size:.9rem}} td,th{{border:1px solid #CFDAD6;padding:.3rem .6rem;text-align:left}}
th{{background:#EAF1EE}} code{{font-family:Consolas,monospace}} .note{{background:#F4F7F6;border-left:4px solid #0F766E;padding:.6rem 1rem}}
</style></head><body>
<h1>PathoPanel report</h1><p>{meta['assay']}. Generated {meta['time']}.</p>{summary}
<h2>Assay</h2><table><tr><th>Oligo</th><th>Sequence (5' to 3')</th><th>Length</th></tr>{oligo_rows}</table>
<h2>How the answer depends on the rule</h2>{rules.to_html() if rules is not None else ""}
<h2>Scorecard by serotype</h2>{card.to_html()}
<h2>Result by serotype</h2>{img("status", "Result by serotype")}
<h2>Mismatch landscape</h2>{img("landscape", "Mismatch landscape")}
<h2>Amplicon length</h2>{img("amplicon", "Amplicon length")}
<h2>Possible redesign</h2>{repair.to_html(index=False)}
<p>Suggestions come from the sequences in the file and need laboratory testing before any use.</p>
<h2>Method and run record</h2><ul>
<li>Each oligo was aligned end to end against both strands (Biopython PairwiseAligner, no gaps inside the oligo); the placement with the fewest mismatches was kept. Binding sites with more than {meta['max_unreadable']} unreadable bases were not called.</li>
<li>Risk score: weight by distance from the 3' end (x4 last base, x3 bases 2 to 3, x2 bases 4 to 5, x1 otherwise) times a pair-type weight (mild 0.5, severe 2.0, other 1.0; pair types after Stadhouders et al., J Mol Diagn 2010;12(1):109-117). The labels and the four-mismatch cutoff are heuristic and not validated against laboratory data.</li>
<li>Serotype was read from FASTA header text. Check it against your metadata.</li>
<li>Reference file: {meta['file']} (SHA-256 {meta['sha']})</li>
<li>Sequences in file: {meta['n_in_file']:,}; removed as shorter than {meta['min_length']:,} bases: {meta['n_short']:,}; analysed: {len(df):,}{sample_note}</li>
<li>Settings: maximum mismatches for a site {meta['max_mm']}; maximum amplicon {meta['max_amp']} bp</li>
<li>Software: Biopython {meta['biopython']}, NumPy {np.__version__}, pandas {pd.__version__}</li>
</ul></body></html>"""
