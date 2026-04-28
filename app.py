import streamlit as st
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
from io import BytesIO

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.pdfgen import canvas

# ------------------ CONFIG ------------------
Entrez.email = "smitangshudas23@gmail.com"

st.set_page_config(page_title="DNA Analyzer", layout="wide")

# ------------------ SESSION INIT ------------------
for key in ["dna", "gc", "counts", "rna", "protein", "gene_results", "blast"]:
    if key not in st.session_state:
        st.session_state[key] = None

# ------------------ DNA INPUT ------------------
st.title("🧬 DNA Analyzer")

dna_input = st.text_area("Enter DNA Sequence")

# ------------------ ANALYZE ------------------
if st.button("Analyze"):
    dna = dna_input.upper()

    if dna:
        seq = Seq(dna)

        counts = {
            "A": dna.count("A"),
            "T": dna.count("T"),
            "G": dna.count("G"),
            "C": dna.count("C")
        }

        gc = ((counts["G"] + counts["C"]) / len(dna)) * 100
        rna = str(seq.transcribe())
        protein = str(seq.translate())

        # SAVE STATE
        st.session_state.dna = dna
        st.session_state.gc = gc
        st.session_state.counts = counts
        st.session_state.rna = rna
        st.session_state.protein = protein

# ------------------ SHOW DNA RESULT ------------------
if st.session_state.dna:
    st.subheader("📊 DNA Analysis")
    st.write(f"Length: {len(st.session_state.dna)}")
    st.write(f"GC Content: {st.session_state.gc:.2f}%")
    st.write(f"Counts: {st.session_state.counts}")
    st.write(f"RNA: {st.session_state.rna}")
    st.write(f"Protein: {st.session_state.protein}")

    import matplotlib.pyplot as plt

    fig1, ax1= plt.subplots()
    ax1.bar(["GC Content"], [st.session_state.gc])
    ax1.set_ylabel("Percentage")
    ax1.set_title("GC Content")
    st.pyplot(fig1)

    counts= st.session_state.counts

    fig2, ax2= plt.subplots()
    ax2.bar(counts.keys(), counts.values())
    ax2.set_title("Nucleotide Distribution")
    st.pyplot(fig2)

# ------------------ NCBI GENE SEARCH ------------------
st.subheader("🧬 NCBI Gene Search")

gene_query = st.text_input("Enter gene name (e.g. BRCA1)")

if st.button("Search Gene"):
    if gene_query:
        with st.spinner("Searching NCBI..."):
            try:
                handle = Entrez.esearch(db="nucleotide", term=gene_query, retmax=3)
                record = Entrez.read(handle)
                handle.close()

                ids = record["IdList"]

                if ids:
                    fetch = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text")
                    data = fetch.read()
                    fetch.close()

                    st.session_state.gene_results = data
                else:
                    st.warning("No results found")

            except Exception as e:
                st.error(f"Error: {e}")

# SHOW GENE RESULTS
if st.session_state.gene_results:
    st.subheader("📄 Gene Results")
    st.text_area("Preview", st.session_state.gene_results[:2000], height=300)

# ------------------ BLAST ------------------
st.subheader("🔬 BLAST")

if st.button("Run BLAST"):
    if st.session_state.dna:
        with st.spinner("Running BLAST..."):
            try:
                result = NCBIWWW.qblast("blastn", "nt", st.session_state.dna)

                blast_records = list(NCBIXML.parse(result))

               if blast_records:
                      st.session_state.blast = blast_records[0]
               else:
                 st.session_state.blast= None
                 st.warning("No BLAST results found")

# SHOW BLAST
if st.session_state.blast:
    st.subheader("🔎 BLAST Results")

    for align in st.session_state.blast.alignments[:2]:
        st.write(align.title)

        for hsp in align.hsps[:1]:
            st.write(f"E-value: {hsp.expect}")

# ------------------ PDF ------------------

def add_watermark(c):
    c.saveState()

    c.setFillGray(0.92)
    c.setFont("Helvetica-Bold", 60)

    width, height = letter

    # Single clean diagonal watermark
            c.translate(width/2, height/2)
            c.rotate(45)

            c.drawCentredString(0, 0, "Smitangshu Bio Lab")
            c.restoreState()


def add_border(c):
    c.setLineWidth(2)
    c.rect(20, 20, 550, 750)

def add_page_design(canvas, doc):
    add_border(canvas)
    add_watermark(canvas)

def make_pdf():
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)
    styles = getSampleStyleSheet()
    story = []

    # DNA
    dna = st.session_state.get("dna", "")

    story.append(Paragraph("DNA Analyzer Report", styles["Heading1"]))
    story.append(Spacer(1, 10))

    if dna:
        story.append(Paragraph(f"Length: {len(dna)}", styles["Normal"]))
        story.append(Paragraph(f"GC: {st.session_state.gc:.2f}%", styles["Normal"]))
        story.append(Paragraph(f"Counts: {st.session_state.counts}", styles["Normal"]))
        story.append(Paragraph(f"RNA: {st.session_state.rna}", styles["Normal"]))
        story.append(Paragraph(f"Protein: {st.session_state.protein}", styles["Normal"]))
    else:
        story.append(Paragraph("No DNA data", styles["Normal"]))

    story.append(Spacer(1, 15))

    # NCBI
    if st.session_state.gene_results:
        story.append(Paragraph("NCBI Gene Results", styles["Heading2"]))
        story.append(Paragraph(st.session_state.gene_results[:1500], styles["Normal"]))

    story.append(Spacer(1, 15))

    # BLAST
    if st.session_state.blast:
        story.append(Paragraph("BLAST Results", styles["Heading2"]))

        for align in st.session_state.blast.alignments[:2]:
            story.append(Paragraph(align.title, styles["Normal"]))
            for hsp in align.hsps [:1]:
                story.append(Paragraph(f"E-value: {hsp.expect}", styles["Normal"]))

    doc.build(story, onFirstPage=add_page_design, onLaterPages=add_page_design)

    buffer.seek(0)
    return buffer

# ------------------ DOWNLOAD ------------------
if st.button("Download Report"):
    pdf = make_pdf()
    st.download_button("Download PDF", pdf, "dna_report.pdf")

# ------------------ FOOTER ------------------
st.markdown("---")
st.markdown("Made with 💕 by Smitangshu")