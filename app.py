import datetime as dt
import hashlib
from pathlib import Path

import Bio
import pandas as pd
import streamlit as st

import core

st.set_page_config(page_title="PathoPanel", page_icon="🧬", layout="wide")

HERE = Path(__file__).parent
DEMO = next((p for p in (HERE / "demo_dengue.fasta", HERE / "data" / "demo_dengue.fasta") if p.exists()), None)
D1 = "TCAATATGCTGAAACGCGCGAGAAACCG"      # Lanciotti et al. 1992, D1, genome position 134-161
D2 = "TTGCACCAACAGTCAATGTCTTCAGGTTC"     # Lanciotti et al. 1992, D2, genome position 616-644

INK, TEAL, MUTED = "#14232B", "#0F766E", "#5C6F76"
TILE = {"A": ("#D3EBDC", "#17603A"), "C": ("#D2E2F6", "#1F4E8C"),
        "G": ("#F7E5BA", "#7A5200"), "T": ("#F5D3D3", "#8E2B2B")}

st.markdown(f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:wght@600;700&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@500&display=swap');
html, body, [class*="st-"]:not([data-testid="stIconMaterial"]):not(.material-icons):not([class*="icon"]) {{ font-family: 'IBM Plex Sans', sans-serif; }}
h1, h2, h3 {{ font-family: 'Bricolage Grotesque', sans-serif !important; color: {INK}; letter-spacing: -0.01em; }}
#MainMenu, footer {{ visibility: hidden; }}
.block-container {{ padding-top: 4.5rem; max-width: 1200px; }}
.stButton > button[kind="primary"] {{ background: {TEAL}; border: none; font-weight: 600; }}
.stButton > button[kind="primary"]:hover {{ background: #0B5F58; }}
.pp-title {{ font-family:'Bricolage Grotesque',sans-serif; font-size:3rem; font-weight:700; color:{INK}; line-height:1; margin:0; }}
.pp-sub {{ color:{MUTED}; font-size:1.05rem; margin:.5rem 0 1.2rem 0; max-width:44rem; }}
.pp-row {{ display:flex; flex-wrap:wrap; align-items:center; gap:.9rem; margin:.45rem 0; }}
.pp-label {{ width:8.5rem; min-width:6rem; font-size:.85rem; color:{MUTED}; flex-shrink:0; }}
.pp-tiles {{ display:flex; flex-wrap:wrap; gap:2px; font-family:'IBM Plex Mono',monospace; }}
.pp-tile {{ width:1.35rem; height:1.7rem; display:flex; align-items:center; justify-content:center; border-radius:3px; font-size:.85rem; }}
.pp-end {{ font-size:.75rem; color:{MUTED}; margin:0 .35rem; }}
.pp-band {{ display:flex; flex-wrap:wrap; border-top:1px solid #CFDAD6; border-bottom:1px solid #CFDAD6; margin:1.2rem 0 1.5rem 0; }}
.pp-stat {{ flex:1; min-width:11rem; padding:1rem 1.2rem; border-left:1px solid #CFDAD6; }}
.pp-stat:first-child {{ border-left:none; padding-left:0; }}
.pp-num {{ font-family:'Bricolage Grotesque',sans-serif; font-size:2.2rem; font-weight:700; line-height:1.1; }}
.pp-cap {{ color:{MUTED}; font-size:.9rem; }}
.pp-headline {{ font-family:'Bricolage Grotesque',sans-serif; font-size:1.5rem; font-weight:600; color:{INK}; line-height:1.35; max-width:50rem; }}
</style>
""", unsafe_allow_html=True)


def tiles(seq: str, label: str) -> str:
    out = []
    for i, ch in enumerate(seq):
        bg, fg = TILE.get(ch, ("#E4E8EA", "#44525A"))
        edge = f"box-shadow: inset 0 -3px 0 {INK};" if len(seq) - i <= 5 else ""
        out.append(f'<span class="pp-tile" style="background:{bg};color:{fg};{edge}">{ch}</span>')
    return (f'<div class="pp-row"><div class="pp-label">{label}</div><span class="pp-end">5\'</span>'
            f'<div class="pp-tiles">{"".join(out)}</div><span class="pp-end">3\'</span></div>')


st.markdown('<p class="pp-title">PathoPanel</p>'
            '<p class="pp-sub">Check how well a diagnostic primer and probe set fits the viral sequences '
            'in a reference collection, and see exactly where and in which serotypes the mismatches fall.</p>',
            unsafe_allow_html=True)

# -------------------------------------------------------------------- sidebar
sb = st.sidebar
sb.header("Assay")
assay_name = sb.text_input("Assay name", "Lanciotti 1992 D1/D2")
fwd = core.clean_oligo(sb.text_input("Forward primer (5' to 3')", D1))
rev = core.clean_oligo(sb.text_input("Reverse primer (5' to 3')", D2))
prb = core.clean_oligo(sb.text_input("Probe (5' to 3', optional)"))
sb.caption("Only A, C, G and T are accepted. Copy sequences from the published paper.")

sb.header("Reference sequences")
upload = sb.file_uploader("FASTA file", type=["fasta", "fa", "fas", "fna", "txt"])
use_demo = sb.checkbox("Use the included demo file", value=DEMO is not None, disabled=DEMO is None)

with sb.expander("Settings"):
    min_length = st.number_input("Skip sequences shorter than (bases)", 0, 20000, 9000, step=500)
    sample_size = st.number_input("Random sample size (0 = use all)", 0, 50000, 0, step=100)
    balance = st.checkbox("Equal numbers per serotype in the sample", value=False, disabled=sample_size == 0)
    max_mm = st.slider("Most mismatches that still count as a binding site", 2, 10, 6)
    max_unreadable = st.slider("Most unreadable bases allowed in a site", 0, 5, 2)
    max_amp = st.number_input("Longest accepted amplicon (bp)", 100, 5000, 2000, step=100)
    expected_amp = st.number_input("Expected amplicon (bp, 0 = skip)", 0, 5000, 511, step=1)
run = sb.button("Run analysis", type="primary", width="stretch")

oligos = {k: v for k, v in {"Fwd": fwd, "Rev": rev, "Probe": prb}.items() if v}

# ------------------------------------------------------------------ run step
if run:
    errors = []
    if not oligos:
        errors.append("Enter at least one primer or probe.")
    for k, v in oligos.items():
        if not core.valid_oligo(v):
            errors.append(f"{core.LABEL[k]} contains letters other than A, C, G and T.")
    if upload is not None:
        raw, fname = upload.getvalue(), upload.name
    elif use_demo and DEMO is not None:
        raw, fname = DEMO.read_bytes(), DEMO.name
    else:
        raw, fname = None, None
        errors.append("Upload a FASTA file or tick the demo file option.")
    for e in errors:
        st.error(e)
    if not errors:
        records = core.read_fasta(raw)
        if not records:
            st.error("No sequences were read. Check that the file is in FASTA format.")
        else:
            records, n_in, n_short = core.filter_and_sample(records, int(min_length), int(sample_size), balance)
            if not records:
                st.error("No sequences are left after the length filter. Lower the minimum length in Settings.")
            else:
                bar = st.progress(0.0, text=f"Comparing {len(oligos)} oligo(s) with {len(records):,} sequences")
                df, site_log = core.analyse(records, oligos, max_mm, int(max_amp), max_unreadable,
                                            progress=lambda f: bar.progress(f, text="Analysing sequences"))
                bar.empty()
                st.session_state["res"] = {
                    "df": df, "site_log": site_log, "oligos": oligos,
                    "figs": {"status": core.fig_status(df, oligos),
                             "landscape": core.fig_landscape(df, site_log, oligos),
                             "amplicon": core.fig_amplicon(df, int(expected_amp))},
                    "meta": {"assay": assay_name, "file": fname, "sha": hashlib.sha256(raw).hexdigest(),
                             "n_in_file": n_in, "n_short": n_short, "min_length": int(min_length),
                             "sample_size": int(sample_size), "max_mm": max_mm, "max_amp": int(max_amp),
                             "max_unreadable": max_unreadable, "expected_amp": int(expected_amp),
                             "time": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
                             "biopython": Bio.__version__}}

res = st.session_state.get("res")

# ---------------------------------------------------------------- empty state
if res is None:
    preview = {core.LABEL[k]: v for k, v in oligos.items() if core.valid_oligo(v)}
    if preview:
        st.markdown("".join(tiles(v, k) for k, v in preview.items()), unsafe_allow_html=True)
        st.caption("The dark bar under the last five bases marks the 3' end, where mismatches matter most.")
    st.info("Choose the demo file or upload a FASTA file in the sidebar, then click Run analysis. "
            "Large files (thousands of genomes) can take a few minutes.")
    st.stop()

# ------------------------------------------------------------------- results
df, site_log, oligos = res["df"], res["site_log"], res["oligos"]
meta, figs = res["meta"], res["figs"]
n = len(df)
rules = core.rules_table(df, oligos)
card = core.scorecard(df, oligos)

st.markdown("".join(tiles(v, core.LABEL[k]) for k, v in oligos.items()), unsafe_allow_html=True)

if rules is not None:
    A, B, C, D = (rules["All sequences"].iloc[i] for i in range(4))
    st.markdown(
        f'<div class="pp-headline">{A}% of {n:,} sequences carry a binding site for every primer, and '
        f'{B}% also match in the last five bases at each 3\' end. Depending on how many mismatches are '
        f'tolerated, compatibility ranges from {D}% to {A}%.</div>', unsafe_allow_html=True)
    st.markdown(f"""
    <div class="pp-band">
      <div class="pp-stat"><div class="pp-num">{n:,}</div><div class="pp-cap">sequences analysed</div></div>
      <div class="pp-stat"><div class="pp-num">{A}%</div><div class="pp-cap">binding site for every primer</div></div>
      <div class="pp-stat"><div class="pp-num" style="color:{TEAL}">{B}%</div><div class="pp-cap">and 3' ends match</div></div>
      <div class="pp-stat"><div class="pp-num" style="color:#C2373B">{D}%</div><div class="pp-cap">strict: 3' ends match and at most 2 mismatches</div></div>
    </div>""", unsafe_allow_html=True)
else:
    likely = core.pct(df["Status"].isin(["Exact", "Minor"]))
    st.markdown(f'<div class="pp-headline">{likely}% of {n:,} sequences have an exact or minor-mismatch '
                f'site for the probe.</div>', unsafe_allow_html=True)

t1, t2, t3, t4, t5 = st.tabs(["Overview", "Mismatch landscape", "3' end", "Variants and redesign", "Data and downloads"])

with t1:
    if rules is not None:
        st.subheader("How much does the answer depend on the rule?")
        st.caption("The same sequences judged under four rules, as the percentage that pass. "
                   "No rule is the single right one, so report the range.")
        st.dataframe(rules, width="stretch")
        st.subheader("Try your own tolerance")
        st.caption("Choose how strict to be. The chart updates at once because it reuses the stored results.")
        c1, c2 = st.columns(2)
        window = c1.slider("Bases at each primer's 3' end that must match", 0, 8, 5)
        limit = c2.slider("Most mismatches allowed in each primer", 0, 10, 4)
        ex = core.explore(df, oligos, window, limit)
        st.write(f"**{ex['All sequences']}%** of sequences pass when the last {window} bases at each 3' end must "
                 f"match and each primer has at most {limit} mismatches.")
        st.bar_chart(ex.drop("All sequences").rename("Share passing (%)"), height=240)
    st.subheader("Scorecard by serotype")
    st.caption("Exact or minor uses simple labels based on mismatch count, position and pair type. "
               "These labels are heuristic and have not been validated against laboratory data.")
    st.dataframe(card, width="stretch")
    st.subheader("Result by serotype")
    st.pyplot(figs["status"])
    st.download_button("Download this figure (PNG)", core.fig_png(figs["status"]), "pathopanel_status.png",
                       "image/png", key="dl_status")

with t2:
    st.subheader("Where mismatches fall")
    st.caption("Each tile is one base of the oligo, with the 5' end on the left. Colour shows the share of "
               "sequences with a binding site that differ at that base. The box marks the last five bases at the 3' end.")
    st.pyplot(figs["landscape"])
    st.download_button("Download this figure (PNG)", core.fig_png(figs["landscape"]), "pathopanel_landscape.png",
                       "image/png", key="dl_landscape")
    if figs["amplicon"] is not None:
        st.subheader("Amplicon length")
        lengths = df["Amplicon length"].dropna()
        msg = f"Median {lengths.median():.0f} bp, range {lengths.min():.0f} to {lengths.max():.0f} bp."
        if meta["expected_amp"]:
            off = int(((lengths - meta["expected_amp"]).abs() > 10).sum())
            msg += f" {off:,} of {len(lengths):,} differ from the expected {meta['expected_amp']} bp by more than 10 bp."
        st.caption(msg)
        st.pyplot(figs["amplicon"])

with t3:
    tp = core.three_prime_summary(df, oligos)
    if not tp:
        st.info("This tab needs at least one primer.")
    for k, s in tp.items():
        st.subheader(core.LABEL[k])
        st.write(f"{s['n']:,} of {n:,} sequences ({100 * s['n'] / n:.1f}%) have a mismatch in the last five bases.")
        if s["n"]:
            c1, c2 = st.columns(2)
            c1.caption("By serotype")
            c1.dataframe(s["by"], width="stretch")
            c2.caption("Most common mismatch positions")
            c2.dataframe(s["positions"], width="stretch")
            st.caption("Most common mismatch patterns. In 27:T/G, position 27 of the primer is T where the template base is G.")
            st.dataframe(s["types"], width="stretch")

with t4:
    thr = st.slider("Add a degenerate base when another base appears in at least this share of sites (%)",
                    1, 30, 5) / 100
    rep = core.redesign(oligos, site_log, thr)
    st.subheader("Possible redesign")
    st.caption("This is a sequence-level idea only. Degenerate bases can change melting temperature, specificity "
               "and primer-dimer risk, so any redesign needs laboratory testing. Very high degeneracy is not practical.")
    st.dataframe(rep, width="stretch", hide_index=True)
    st.subheader("Alternative: a small mix of serotype-specific versions")
    st.caption("Instead of one highly degenerate oligo, use the most common version from each serotype. "
               "A mix adds cost and can raise primer-dimer risk, so it also needs laboratory testing.")
    for k, o in oligos.items():
        mix = core.serotype_mix(k, o, site_log)
        if mix is None:
            continue
        table, covered, now = mix
        st.write(f"**{core.LABEL[k]}:** a mix of {len(table)} version(s) would match exactly {covered}% of sites, "
                 f"against {now}% for the current oligo.")
        st.dataframe(table, width="stretch", hide_index=True)
    st.subheader("Most common binding-site variants")
    pick = st.selectbox("Oligo", list(oligos), format_func=lambda k: core.LABEL[k])
    vt = core.variant_table(site_log, pick)
    if vt is None:
        st.warning("No binding sites were found for this oligo.")
    else:
        st.caption("Sites are written in the orientation of the oligo, so a perfect match reads like the oligo itself.")
        st.dataframe(vt, width="stretch")

with t5:
    st.subheader("Results table")
    shown = st.multiselect("Show these whole-assay statuses", core.STATUS_ORDER, default=core.STATUS_ORDER)
    st.dataframe(df[df["Status"].isin(shown)], width="stretch", height=380, hide_index=True)
    report_figs = {k: core.fig_b64(f) for k, f in figs.items() if f is not None}
    report = core.build_report(meta, oligos, df, rules if rules is not None else None, card, rep, report_figs) \
        if rules is not None else None
    c1, c2 = st.columns(2)
    c1.download_button("Download results table (CSV)", df.to_csv(index=False).encode(),
                       "pathopanel_results.csv", "text/csv", width="stretch")
    if report:
        c2.download_button("Download report (HTML)", report.encode("utf-8"),
                           "pathopanel_report.html", "text/html", width="stretch")
    st.caption(f"Run record: {meta['file']} (SHA-256 {meta['sha'][:16]}...), {meta['n_in_file']:,} sequences in file, "
               f"{meta['n_short']:,} removed as shorter than {meta['min_length']:,} bases, {n:,} analysed, "
               f"run {meta['time']}.")

with st.expander("Method, limits and references"):
    st.markdown("""
**Method.** Each oligo is aligned end to end against both strands of every sequence, with no gaps inside the oligo,
and the placement with the fewest mismatches is kept. Serotype is read from the FASTA header text.

**What this shows.** Sequence compatibility between oligos and reference sequences. A match does not prove
amplification, and a mismatch does not prove diagnostic failure.

**Limits.** The mismatch labels, the four-mismatch cutoff and the weights are heuristic and have not been validated
against laboratory data. The cutoff follows a synthetic-template study in which four mismatches in a 20-base primer blocked
the reaction; here it is applied to any primer and any position. Pair-type weights come from a study of 3' end mismatches
and are applied here to all positions.
Public sequence collections are not a random sample of circulating virus, and repeated sequences from single
studies can inflate some groups.

**References.**
- Lanciotti RS, et al. 1992. Journal of Clinical Microbiology (D1 and D2 primers).
- Lefever S, Pattyn F, Hellemans J, Vandesompele J. 2013. Single-nucleotide polymorphisms and other mismatches reduce
  performance of quantitative PCR assays. Clinical Chemistry 59:1470-1480.
- Vanneste K, Garlant L, Broeders S, et al. 2018. Application of whole genome data for in silico evaluation of primers and
  probes routinely employed for the detection of viral species by RT-qPCR using dengue virus as a case study. BMC
  Bioinformatics 19:312.
- Otwell T, et al. 2025. Reality check: testing the in silico predictions of false negative results due to mutations in
  SARS-CoV-2 PCR assays using templates with mismatches in vitro. Frontiers in Cellular and Infection Microbiology 15:1524025.
- Stadhouders R, et al. 2010. The effect of primer-template mismatches on the detection and quantification of nucleic acids
  using the 5' nuclease assay. Journal of Molecular Diagnostics 12(1):109-117. doi:10.2353/jmoldx.2010.090035
""")
