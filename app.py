import streamlit as st
from Bio.Seq import Seq
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from io import BytesIO

Entrez.email = "smitangshudas23@gmail.com"

st.set_page_config(page_title="DNA Analyzer", layout="wide")

# ---------------- STATE ----------------
if "analysis" not in st.session_state:
    st.session_state.analysis = None
if "gene" not in st.session_state:
    st.session_state.gene = None
if "blast" not in st.session_state:
    st.session_state.blast = None

# ---------------- INPUT ----------------
st.title("🧬 DNA Analyzer")

dna_input = st.text_area("Enter DNA Sequence")

# CLEAN INPUT
def clean_dna(seq):
    return seq.replace("\n", "").replace(" ", "").upper()

def valid_dna(seq):
    return all(c in "ATGC" for c in seq)

# ---------------- ANALYZE ----------------
if st.button("Analyze"):
    dna = clean_dna(dna_input)

    if not dna or not valid_dna(dna):
        st.error("Invalid DNA sequence")
    else:
        seq = Seq(dna)

        length = len(dna)
        g = dna.count("G")
        c = dna.count("C")
        gc = round(((g+c)/length)*100, 2)

        counts = {
            "A": dna.count("A"),
            "T": dna.count("T"),
            "G": dna.count("G"),
            "C": dna.count("C"),
        }

        rna = seq.transcribe()
        protein = seq.translate()

        st.session_state.analysis = {
            "dna": dna,
            "length": length,
            "gc": gc,
            "counts": counts,
            "rna": rna,
            "protein": protein
        }

# ---------------- SHOW ANALYSIS ----------------
if st.session_state.analysis:
    data = st.session_state.analysis

    st.subheader("📊 DNA Analysis")
    st.write(f"Length: {data['length']}")
    st.write(f"GC Content: {data['gc']}%")
    st.write(data["counts"])
    st.write("RNA:", data["rna"])
    st.write("Protein:", data["protein"])

    fig, ax = plt.subplots()
    ax.bar(data["counts"].keys(), data["counts"].values())
    st.pyplot(fig)

# ---------------- NCBI GENE SEARCH ----------------
st.subheader("🧬 NCBI Gene Search")

gene_query = st.text_input("Enter gene name (e.g. BRCA1)")

if st.button("Search Gene"):
    if gene_query:
        handle = Entrez.esearch(db="nucleotide", term=gene_query, retmax=3)
        record = Entrez.read(handle)
        ids = record["IdList"]

        if ids:
            fetch = Entrez.efetch(
                db="nucleotide",
                id=",".join(ids),
                rettype="fasta",
                retmode="text"
            )
            st.session_state.gene = fetch.read()

# SHOW GENE RESULT
if st.session_state.gene:
    st.subheader("📄 Gene Results")
    st.text_area("Preview", st.session_state.gene[:2000], height=300)

# ---------------- BLAST ----------------
st.subheader("🔬 BLAST")

def quality(e):
    if e < 1e-5:
        return "🟢 Strong"
    elif e < 0.01:
        return "🟡 Moderate"
    return "🔴 Weak"

if st.button("Run BLAST"):
    if st.session_state.analysis:
        dna = st.session_state.analysis["dna"][:200]

        try:
            with st.spinner("Running BLAST..."):
                result = NCBIWWW.qblast("blastn", "nt", dna)
                blast_record = NCBIXML.read(result)

                st.session_state.blast = blast_record
                st.success("BLAST completed")

        except Exception as e:
            st.error(f"BLAST failed: {e}")

# SHOW BLAST
if st.session_state.blast:
    st.subheader("🔍 BLAST Results")

    count = 0
    for align in st.session_state.blast.alignments:
        for hsp in align.hsps:

            if hsp.expect > 0.01:
                continue

            count += 1
            if count > 3:
                break

            with st.expander(align.title[:60]):
                st.write(f"Score: {hsp.score}")
                st.write(f"E-value: {hsp.expect}")
                st.write(f"Quality: {quality(hsp.expect)}")

# ---------------- PDF ----------------
def add_watermark(c):
    c.saveState()

    c.setFillGray(0.97)
    c.translate(300, 400)
    c.rotate(45)

    c.drawCentredString(0, 0, "Smitangshu Bio Lab")
    c.drawCentredString(200, 200, "Smitangshu Bio Lab")

    c.restoreState()

def add_border(c):
    c.setLineWidth(2)
    c.rect(30, 30, 530, 730)

def make_pdf():
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)
    styles = getSampleStyleSheet()
    story = []

    story.append(Paragraph("DNA Analyzer Report", styles['Heading1']))
    story.append(Spacer(1, 10))

    if st.session_state.analysis:
        a = st.session_state.analysis
        story.append(Paragraph(f"Length: {a['length']}", styles['Normal']))
        story.append(Paragraph(f"GC: {a['gc']}%", styles['Normal']))
        story.append(Paragraph(str(a["counts"]), styles['Normal']))
        story.append(Paragraph(str(a["rna"]), styles['Normal']))
        story.append(Paragraph(str(a["protein"]), styles['Normal']))

    if st.session_state.gene:
        story.append(Paragraph("Gene Results", styles['Heading2']))
        story.append(Paragraph(st.session_state.gene[:1000], styles['Normal']))

    if st.session_state.blast:
        story.append(Paragraph("BLAST Results", styles['Heading2']))
        for align in st.session_state.blast.alignments[:2]:
            for hsp in align.hsps:
                story.append(Paragraph(align.title, styles['Normal']))
                story.append(Paragraph(f"E-value: {hsp.expect}", styles['Normal']))

    def add_page_design(canvas, doc):
        add_border(canvas)
        add_watermark(canvas)

    doc.build(
        story,
        onFirstPage=add_page_design,
        onLaterPages=add_page_design)

    buffer.seek(0)
    return buffer

if st.button("Download Report"):
    pdf = make_pdf()
    st.download_button("Download PDF", pdf, "dna_report.pdf")

# ---------------- FOOTER ----------------
st.markdown("---")
st.markdown("Made with 💕 by Smitangshu")