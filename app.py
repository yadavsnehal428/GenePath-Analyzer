import pandas as pd
import streamlit as st

from genepath import __version__
from genepath.pipeline import annotate, clean_dna
from genepath.reference import accession_from_header, fetch_transcript, find_transcript_by_blast
from genepath.vep import annotate_rows

st.set_page_config(page_title="GenePath Analyzer", page_icon="🧬", layout="wide")


@st.cache_data(show_spinner=False, ttl=24 * 3600)
def get_transcript(accession: str):
    return fetch_transcript(accession)


@st.cache_data(show_spinner=False, ttl=24 * 3600)
def get_vep(rows: tuple, accession: str):
    return annotate_rows([dict(r) for r in rows], accession)


# ---------------- Sidebar ----------------
st.sidebar.title("🧬 GenePath Analyzer")
st.sidebar.markdown(
    "Finds the variants in a human transcript sequence and describes them in "
    "**HGVS notation** with their effect on the protein."
)
st.sidebar.markdown("""
**How to use**
1. Upload or paste a DNA/cDNA sequence (FASTA or plain text)
2. Give the RefSeq transcript it belongs to (e.g. `NM_007194.4`). It is read from the FASTA header when present.
3. Click **Analyse**
""")
st.sidebar.markdown("---")
st.sidebar.caption(f"Version {__version__}")
st.sidebar.caption("Mentored by Dr. Kushagra Kashyap")
st.sidebar.caption("Developed by Snehal Yadav")
st.sidebar.caption("For research and education only, not for clinical use.")

# ---------------- Input ----------------
st.title("🧬 GenePath Analyzer")
st.caption("Variant detection and HGVS annotation against RefSeq transcripts")

col1, col2 = st.columns(2)
with col1:
    uploaded = st.file_uploader("Upload sequence (.fasta, .fa, .txt)", type=["fasta", "fa", "txt"])
with col2:
    pasted = st.text_area("…or paste the sequence", height=120)

raw = uploaded.read().decode("utf-8") if uploaded else pasted
header = raw.strip().splitlines()[0] if raw.strip().startswith(">") else ""

accession = st.text_input(
    "Reference transcript (RefSeq accession)",
    value=accession_from_header(header) or "",
    placeholder="NM_007194.4",
    help="Leave empty to identify the transcript by BLAST (slower).",
)
use_vep = st.checkbox("Annotate with Ensembl VEP (SIFT, PolyPhen, ClinVar, gnomAD)", value=True)

if st.button("🚀 Analyse", type="primary"):
    if not raw.strip():
        st.warning("Please provide a sequence.")
        st.stop()
    try:
        query = clean_dna(raw)
    except ValueError as e:
        st.error(str(e))
        st.stop()

    if not accession:
        with st.spinner("No accession given: searching human RefSeq mRNAs with BLAST…"):
            accession = find_transcript_by_blast(query)
        if not accession:
            st.error("BLAST found no matching human RefSeq transcript.")
            st.stop()
        st.info(f"Best matching transcript: **{accession}**")

    with st.spinner(f"Fetching {accession} from NCBI…"):
        try:
            tx = get_transcript(accession)
        except Exception as e:
            st.error(f"Could not fetch {accession}: {e}")
            st.stop()

    with st.spinner("Aligning and annotating…"):
        rows = annotate(tx, query)
    if use_vep and rows:
        with st.spinner("Querying Ensembl VEP…"):
            try:
                rows = get_vep(tuple(tuple(r.items()) for r in rows), tx.accession)
            except Exception as e:
                st.warning(f"Ensembl VEP could not be reached ({e}). Showing HGVS results only.")
    df = pd.DataFrame(rows)
    if "VEP note" in df.columns and df["VEP note"].astype(bool).any():
        notes = df.loc[df["VEP note"].astype(bool), "VEP note"]
        st.warning(f"VEP could not annotate {len(notes)} of {len(df)} variants. "
                   f"First error: {notes.iloc[0]}")

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Gene", tx.gene)
    m2.metric("Transcript", tx.accession)
    m3.metric("Coding region", f"{tx.cds_start}–{tx.cds_end}")
    m4.metric("Variants found", len(df))

    if df.empty:
        st.success("No differences from the reference transcript in the region covered.")
        st.stop()

    tab1, tab2 = st.tabs(["Variant table", "Summary"])
    with tab1:
        st.dataframe(
            df, use_container_width=True, hide_index=True,
            column_config={"gnomAD AF": st.column_config.NumberColumn(format="%.2e")},
        )
        if "SIFT" in df.columns:
            st.caption(
                "SIFT (≤0.05 = deleterious) and PolyPhen-2 (closer to 1 = damaging) are "
                "shown for missense variants only. ClinVar significance is for the same "
                "allele where Ensembl reports it. An empty gnomAD AF means the variant "
                "was not seen in gnomAD, which is typical for rare or novel variants."
            )
        st.download_button("💾 Download CSV", df.to_csv(index=False),
                           file_name=f"GenePath_{tx.gene}_{tx.accession}.csv")
    with tab2:
        st.subheader("Consequences")
        st.bar_chart(df["Consequence"].value_counts())
        if "Impact" in df.columns and df["Impact"].astype(bool).any():
            st.subheader("VEP impact")
            st.bar_chart(df["Impact"].replace("", "not annotated").value_counts())
        st.caption(
            "Consequence terms follow the Sequence Ontology. Protein changes are "
            "predicted from the transcript sequence (p. in parentheses, as HGVS "
            "recommends for predicted changes)."
        )
