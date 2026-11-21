# PathoPanel

Checks how well a diagnostic primer and probe set fits the sequences in a FASTA file, split by dengue serotype.
It reports binding-site coverage, 3' end mismatches, a mismatch landscape, common site variants and a redesign idea.

The app measures sequence compatibility only. A match does not prove amplification, and a mismatch does not prove
diagnostic failure.

Run locally: `pip install -r requirements.txt` then `streamlit run app.py`.
Put a small FASTA file named `demo_dengue.fasta` next to `app.py` to enable the demo option.

Primers D1 and D2 are from Lanciotti et al. 1992, Journal of Clinical Microbiology.
