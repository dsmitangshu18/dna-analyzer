import streamlit as st
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
import matplotlib.pyplot as plt
from io import BytesIO

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.pdfgen import canvas

# ================= CONFIG =================
Entrez.email = "smitangshudas23@gmail.com"
st.set_page_config(page_title="DNA Analyzer Pro", layout="wide")

# ================= SESSION =================
for key in ["dna", "protein", "rna", "blast", "ncbi_results"]:
    if key not in st.session_state:
        st.session_state[key] = None

# ================= INPUT =================
st.title("🧬 DNA Analyzer Pro")

dna_input = st.text_area("Enter DNA Sequence")

# ================= ANALYZE BUTTON =================
if st.button("🔬 Analyze DNA"):
    dna = dna_input.upper().replace("\n", "").strip()

    if dna:
        st.session_state.dna = dna

        seq = Seq(dna)
        st.session_state.rna = str(seq.transcribe())
        st.session_state.protein = str(seq.translate())

        st.success("Analysis Complete!")

# ================= TABS =================
tab1, tab2, tab3, tab4 = st.tabs([
    "🧬 DNA Analysis",
    "🧪 BLAST",
    "🔎 NCBI",
    "📄 Report"
])

# ================= DNA ANALYSIS =================
with tab1:
    if st.session_state.dna:
        dna = st.session_state.dna

        st.subheader("DNA Analysis")

        st.write(f"Length: {len(dna)}")

        # GC
        gc = (dna.count("G") + dna.count("C")) / len(dna) * 100
        st.write(f"GC Content: {gc:.2f}%")

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

        # RNA + Protein
        st.write(f"RNA: {st.session_state.rna}")
        st.write(f"Protein: {st.session_state.protein}")

# ================= BLAST =================
with tab2:
    if st.session_state.dna:

        if st.button("Run BLAST"):
            with st.spinner("Running BLAST..."):
                try:
                    result = NCBIWWW.qblast("blastn", "nt", st.session_state.dna)
                    records = list(NCBIXML.parse(result))

                    st.session_state.blast = records[0]
                    st.success("BLAST Completed")

                except:
                    st.warning("Using Fake BLAST")

                    class FakeHSP:
                        def __init__(self):
                            self.expect = 0.0001
                            self.score = 50
                            self.identities = 20
                            self.align_length = 25

                    class FakeAlign:
                        def __init__(self):
                            self.title = "Fake Gene Match - Homo sapiens"
                            self.hsps = [FakeHSP()]

                    class FakeBlast:
                        def __init__(self):
                            self.alignments = [FakeAlign()]

                    st.session_state.blast = FakeBlast()

        # SHOW RESULT
        if st.session_state.blast:
            for align in st.session_state.blast.alignments[:3]:
                for hsp in align.hsps[:1]:

                    identity = (hsp.identities / hsp.align_length) * 100

                    st.markdown(f"### 🧬 {align.title}")
                    st.write(f"E-value: {hsp.expect}")
                    st.write(f"Score: {hsp.score}")
                    st.write(f"Identity: {identity:.2f}%")
                    st.write(f"Length: {hsp.align_length}")

# ================= NCBI =================
with tab3:
    gene_name = st.text_input("Search gene (e.g. BRCA1)")

    if st.button("Search NCBI"):
        if gene_name:
            with st.spinner("Searching..."):
                handle = Entrez.esearch(db="gene", term=gene_name, retmax=5)
                record = Entrez.read(handle)

                st.session_state.ncbi_results = record["IdList"]

            st.success("Results loaded!")

    if st.session_state.ncbi_results:
        for gene in st.session_state.ncbi_results:
            st.write(f"Gene ID: {gene}")

# ================= PDF =================

def add_watermark(c, doc):
    c.saveState()
    c.setFont("Helvetica-Bold", 30)
    c.setFillGray(0.85)
    c.translate(300, 400)
    c.rotate(45)
    c.drawCentredString(0, 0, "DNA ANALYZER PRO")
    c.restoreState()

def add_border(c, doc):
    c.saveState()
    c.setLineWidth(2)
    c.rect(30, 30, 535, 750)
    c.restoreState()

def add_design(c, doc):
    add_border(c, doc)
    add_watermark(c, doc)

def make_pdf():
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer)

    styles = getSampleStyleSheet()
    story = []

    # DNA
    if st.session_state.dna:
        dna = st.session_state.dna

        story.append(Paragraph("DNA Analysis", styles["Heading2"]))
        story.append(Paragraph(f"Sequence: {dna}", styles["Normal"]))
        story.append(Paragraph(f"Length: {len(dna)}", styles["Normal"]))

        gc = (dna.count("G") + dna.count("C")) / len(dna) * 100
        story.append(Paragraph(f"GC Content: {gc:.2f}%", styles["Normal"]))

        story.append(Paragraph(f"RNA: {st.session_state.rna}", styles["Normal"]))
        story.append(Paragraph(f"Protein: {st.session_state.protein}", styles["Normal"]))

    # NCBI
    if st.session_state.ncbi_results:
        story.append(Spacer(1, 10))
        story.append(Paragraph("NCBI Results", styles["Heading2"]))

        for gene in st.session_state.ncbi_results:
            story.append(Paragraph(f"Gene ID: {gene}", styles["Normal"]))

    # BLAST
    if st.session_state.blast:
        story.append(Spacer(1, 10))
        story.append(Paragraph("BLAST Results", styles["Heading2"]))

        for align in st.session_state.blast.alignments[:3]:
            for hsp in align.hsps[:1]:

                identity = (hsp.identities / hsp.align_length) * 100

                story.append(Paragraph(f"{align.title}", styles["Normal"]))
                story.append(Paragraph(f"E-value: {hsp.expect}", styles["Normal"]))
                story.append(Paragraph(f"Identity: {identity:.2f}%", styles["Normal"]))

    doc.build(story, onFirstPage=add_design, onLaterPages=add_design)

    buffer.seek(0)
    return buffer

# ================= DOWNLOAD =================
with tab4:
    if st.session_state.dna:
        pdf = make_pdf()

        st.download_button(
            label="📄 Download Report",
            data=pdf,
            file_name="dna_report.pdf",
            mime="application/pdf"
        )

# ------------------ FOOTER ------------------
st.markdown("---")
st.markdown("Made with 💕 by Smitangshu")