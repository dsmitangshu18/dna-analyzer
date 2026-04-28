import streamlit as st
from Bio.Seq import Seq
from Bio import Entrez
import pandas as pd
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.pdfgen import canvas
from datetime import datetime
import io

# Email for NCBI
Entrez.email = "smitangshudas23@gmail.com"

st.set_page_config(page_title="DNA Analyzer", layout="wide")

st.title("🧬 DNA Sequence Analyzer")
st.caption("Professional Bioinformatics Tool")

# ---------------- INPUT ----------------
dna_input = st.text_area("Enter DNA Sequence")

# ---------------- SESSION ----------------
if "analysis_done" not in st.session_state:
    st.session_state.analysis_done = False

if "ncbi_results" not in st.session_state:
    st.session_state.ncbi_results = ""

# ---------------- ANALYZE BUTTON ----------------
if st.button("Analyze"):
    if dna_input:
        dna = dna_input.upper()
        if set(dna).issubset({"A","T","G","C"}):
            st.session_state.analysis_done = True
            st.session_state.dna = dna
        else:
            st.error("Invalid DNA sequence")
    else:
        st.warning("Enter DNA first")

# ---------------- ANALYSIS ----------------
if st.session_state.analysis_done:

    dna = st.session_state.dna

    st.subheader("🧪 Sequences")
    st.code(dna)

    seq = Seq(dna)
    reverse_comp = seq.reverse_complement()
    rna = seq.transcribe()
    protein = seq.translate()

    st.code(reverse_comp)
    st.code(rna)
    st.code(protein)

    # GC Content
    gc = (dna.count("G") + dna.count("C")) / len(dna) * 100
    st.success(f"GC Content: {gc:.2f}%")

    # Nucleotide graph
    counts = {
        "A": dna.count("A"),
        "T": dna.count("T"),
        "G": dna.count("G"),
        "C": dna.count("C")
    }

    fig, ax = plt.subplots()
    ax.bar(counts.keys(), counts.values())
    st.pyplot(fig)

# ---------------- NCBI SEARCH ----------------
st.subheader("🌐 NCBI Sequence Search")

ncbi_query = st.text_input("Enter DNA / Gene for NCBI search")

if st.button("Search NCBI Database"):

    if not ncbi_query:
        st.warning("Enter something to search")
    else:
        try:
            with st.spinner("Searching NCBI..."):

                handle = Entrez.esearch(
                    db="nucleotide",
                    term=ncbi_query,
                    retmax=3
                )
                record = Entrez.read(handle)
                handle.close()

                ids = record["IdList"]

                if not ids:
                    st.warning("No matches found")
                else:
                    st.success(f"Found {len(ids)} matches")

                    fetch = Entrez.efetch(
                        db="nucleotide",
                        id=",".join(ids),
                        rettype="fasta",
                        retmode="text"
                    )

                    results = fetch.read()
                    fetch.close()

                    st.session_state.ncbi_results = results

        except Exception as e:
            st.error(f"Error: {e}")

# ---------------- SHOW NCBI RESULTS ----------------
if st.session_state.ncbi_results:
    st.subheader("📄 Retrieved sequences")
    st.text_area("Preview", st.session_state.ncbi_results[:2000], height=300)

# ---------------- PDF GENERATION ----------------
def create_pdf():
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4)
    styles = getSampleStyleSheet()

    content = []

    content.append(Paragraph("DNA Analysis Report", styles["Title"]))
    content.append(Spacer(1, 10))

    if st.session_state.analysis_done:
        dna = st.session_state.dna
        content.append(Paragraph(f"DNA: {dna}", styles["Normal"]))
        content.append(Spacer(1, 10))

    if st.session_state.ncbi_results:
        content.append(Paragraph("NCBI Results:", styles["Heading2"]))
        content.append(Paragraph(st.session_state.ncbi_results[:1000], styles["Normal"]))

    doc.build(content)

    buffer.seek(0)
    return buffer

# ---------------- DOWNLOAD BUTTON ----------------
if st.button("Download Professional Report"):
    pdf = create_pdf()
    st.download_button(
        label="Download PDF",
        data=pdf,
        file_name="dna_report.pdf",
        mime="application/pdf"
    )

# ---------------- FOOTER ----------------
st.markdown("---")
st.markdown("❤️ Made with 💕 by Smitangshu")