import streamlit as st
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import pairwise2
import matplotlib.pyplot as plt
import pandas as pd
import re
from io import StringIO
import requests

# -----------------------------
# PAGE CONFIG
# -----------------------------
st.set_page_config(
    page_title="GenePath Analyzer",
    page_icon="🧬",
    layout="wide"
)

# -----------------------------
# CUSTOM CSS
# -----------------------------
st.markdown("""
<style>
.main-title {font-size:40px; font-weight:bold; color:#2E7D32; margin-bottom:20px;}
.sub-header {font-size:20px; font-weight:bold; color:#1565C0; margin-top:10px;}
.report-box {background:#f9f9f9; padding:15px; border-radius:8px; margin-top:10px;}
</style>
""", unsafe_allow_html=True)

# -----------------------------
# SIDEBAR
# -----------------------------
st.sidebar.title("🧬 GenePath Analyzer")
page = st.sidebar.radio("Navigation", ["🔬 Sequence Analysis", "ℹ️ About"])
st.sidebar.markdown("---")

st.sidebar.markdown("""
### 📌 Instructions
- Upload a DNA sequence file (`.txt` or `.fasta`)
- Or paste sequence manually
- Only valid nucleotides allowed: **A, T, G, C**
- Click **Analyze Sequence**
- Download results as CSV
""")

st.sidebar.markdown("---")
st.sidebar.caption("Version 1.4")
st.sidebar.caption("Developed by Snehal Yadav") 
st.sidebar.caption("Mentored by Dr. Kushagra Kashyap")

# -----------------------------
# FUNCTIONS
# -----------------------------
def extract_gene(text):
    match = re.search(r"\((.*?)\)", text)
    if match:
        return match.group(1)
    for word in text.split():
        if word.isupper() and len(word) <= 10:
            return word
    return "Unknown"

# Expanded local database
local_db = {
    "BRCA1": "DNA Repair Pathway",
    "TP53": "Cell Cycle Regulation",
    "EGFR": "Signal Transduction",
    "KRAS": "MAPK Pathway",
    "MYC": "Cell Growth Regulation",
    "PTEN": "PI3K/AKT Pathway",
    "APC": "Wnt Signaling Pathway",
    "RB1": "Cell Cycle Control",
    "ALK": "Tyrosine Kinase Signaling",
    "BRAF": "MAPK/ERK Pathway",
    "CDKN2A": "Cell Cycle Arrest",
    "JAK2": "JAK/STAT Pathway",
    "NOTCH1": "Notch Signaling Pathway",
    "SMAD4": "TGF-beta Signaling",
    "VHL": "Hypoxia Response Pathway",
    "NF1": "RAS Regulation",
    "PIK3CA": "PI3K/AKT Pathway",
    "FGFR1": "FGF Signaling Pathway",
    "CTNNB1": "Wnt/Beta-Catenin Pathway",
    "MLH1": "Mismatch Repair Pathway",
    "MSH2": "Mismatch Repair Pathway",
    "ATM": "DNA Damage Response",
    "CHEK2": "Cell Cycle Checkpoint",
    "RET": "Tyrosine Kinase Signaling",
    "TERT": "Telomere Maintenance"
}

# KEGG API lookup
def get_pathway_kegg(gene):
    try:
        url = f"http://rest.kegg.jp/link/pathway/hsa:{gene}"
        response = requests.get(url)
        if response.status_code == 200 and response.text.strip():
            pathways = [line.split("\t")[1] for line in response.text.strip().split("\n")]
            return pathways
        else:
            return []
    except Exception:
        return []

def get_pathway_name(pathway_id):
    try:
        url = f"http://rest.kegg.jp/get/{pathway_id}"
        response = requests.get(url)
        if response.status_code == 200:
            for line in response.text.split("\n"):
                if line.startswith("NAME"):
                    return line.replace("NAME", "").strip()
    except Exception:
        pass
    return "Unknown Pathway"

def get_pathway(gene):
    if gene in local_db:
        return [local_db[gene]]
    pathways = get_pathway_kegg(gene)
    if pathways:
        return [get_pathway_name(p) for p in pathways]
    return []

def classify_mutation(ref, query):
    if ref == "-":
        return "Insertion"
    elif query == "-":
        return "Deletion"
    else:
        return "Substitution"

def predict_harmfulness(mutation_df):
    n = len(mutation_df)
    if n == 0:
        return "No mutation detected ✅"
    elif n <= 2:
        return "Likely benign 🟢"
    elif n <= 5:
        return "Possibly harmful 🟠"
    else:
        return "Likely pathogenic 🔴"

def plot_mutation_graph(mutation_df):
    if mutation_df.empty:
        st.write("No mutations to plot")
        return

    counts = mutation_df.groupby('Type').size()
    types = ['Substitution', 'Insertion', 'Deletion']
    values = [counts.get(t, 0) for t in types]
    colors = ['#FF7043', '#42A5F5', '#66BB6A']

    fig, ax = plt.subplots(figsize=(6,3))  # FIXED SIZE
    bars = ax.bar(types, values, color=colors, edgecolor='black')

    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height}', xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0,3), textcoords="offset points", ha='center', va='bottom',
                    fontsize=10, fontweight='bold')

    ax.set_ylabel("Mutations", fontsize=10)
    ax.set_xlabel("Type", fontsize=10)
    ax.set_title("Mutation Distribution", fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    st.pyplot(fig)

# -----------------------------
# ABOUT PAGE
# -----------------------------
if page == "ℹ️ About":
    st.title("ℹ️ About GenePath Analyzer")

    st.markdown("""
        **GenePath Analyzer** is an interactive bioinformatics application designed to analyze DNA sequences, detect mutations, and map their biological significance in a clear and user-friendly manner.

   
    ---
    ### 🚀 Key Features
    - 🔍 Real-time NCBI BLAST integration  
    - 🧾 Detailed mutation table with classification (Insertion, Deletion, Substitution)  
    - 📊 Professional mutation visualization (stacked graph)  
    - 🧬 Gene identification and pathway mapping  
    - ⬇️ Downloadable results in CSV format  

    ---
    ### 🛠️ Technology Stack
    - **Python**
    - **Streamlit**
    - **Biopython**
    - **Pandas**
    - **Matplotlib**

    ---
    ### ⚠️ Disclaimer
    This application is intended for **educational and research purposes only**.  
    It should **not be used for clinical or diagnostic decisions**.

    ---
    ### 👨‍💻 Developer Note
    This tool is designed to simplify complex bioinformatics workflows and make genetic analysis more accessible to students and researchers.
    """)
#------------------------------
# MAIN ANALYZER
# -----------------------------
if page == "🔬 Sequence Analysis":
    st.markdown('<div class="main-title">🧬 GenePath Analyzer</div>', unsafe_allow_html=True)

    uploaded_file = st.file_uploader("Upload DNA sequence (.txt or .fasta)", type=["txt","fasta"])
    sequence = ""
    if uploaded_file is not None:
        sequence = uploaded_file.read().decode("utf-8")
        st.success("File uploaded successfully ✅")

    sequence_input = st.text_area("Paste DNA sequence")
    if sequence_input:
        sequence = sequence_input

    if st.button("🚀 Analyze Sequence"):
        if not sequence:
            st.warning("Please provide a DNA sequence")
        else:
            with st.spinner("Performing BLAST search via NCBI... ⏳"):
                try:
                    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
                    blast_record = NCBIXML.read(result_handle)

                    if not blast_record.alignments:
                        st.error("No matches found")
                        st.stop()

                    alignment = blast_record.alignments[0]
                    hsp = alignment.hsps[0]
                    ref_seq = hsp.sbjct
                    query_seq = hsp.query
                    gene_info = alignment.title

                    alignments = pairwise2.align.globalxx(ref_seq, query_seq)
                    ref_aligned, query_aligned, _, _, _ = alignments[0]

                    mutation_list = []
                    for i in range(len(ref_aligned)):
                        if ref_aligned[i] != query_aligned[i]:
                            mutation_list.append({
                                "Position": i+1,
                                "Reference": ref_aligned[i],
                                "Query": query_aligned[i],
                                "Type": classify_mutation(ref_aligned[i], query_aligned[i])
                            })
                    mutation_df = pd.DataFrame(mutation_list)

                    gene = extract_gene(gene_info)
                    pathways = get_pathway(gene)

                    tabs = st.tabs(["Gene Information","Mutation Table","Mutation Graph"])

                    with tabs[0]:
                        st.subheader("🧬 Gene")
                        st.info(gene)

                        st.subheader("⚠️ Pathogenicity Assessment")
                        st.success(predict_harmfulness(mutation_df))

                        st.subheader("🧪 Pathway")
                        if pathways:
                            for p in pathways:
                                st.markdown(f"This gene is associated with the **{p}**.")
                        else:
                            st.warning("Pathway information not available in current sources.")

                        st.subheader("🔍 BLAST Match")
                        st.code(gene_info[:300])  # CLEAN DISPLAY

                    with tabs[1]:
                        st.subheader("🧬 Mutation Table")
                        if mutation_df.empty:
                            st.write("No mutations detected")
                        else:
                            st.dataframe(mutation_df)

                    with tabs[2]:
                        st.subheader("📊 Mutation Distribution")
                        plot_mutation_graph(mutation_df)

                    csv_buffer = StringIO()
                    mutation_df.to_csv(csv_buffer, index=False)
                    if pathways:
                        csv_buffer.write("\nPathways\n" + "; ".join(pathways) + "\n")

                    st.download_button("💾 Download CSV Report", csv_buffer.getvalue(),
                                       file_name="GenePath_Report.csv")

                except Exception as e:
                    st.error(f"Error: {e}")