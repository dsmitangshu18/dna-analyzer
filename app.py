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
for key in ["dna", "protein", "rna", "blast", "ncbi_results", "orfs"]:
    if key not in st.session_state:
        st.session_state[key] = None

# ================= ORF FUNCTION =================
def find_orfs(dna):
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []

    for frame in range(3):
        i = frame
        while i < len(dna) - 2:
            codon = dna[i:i+3]

            if codon == "ATG":
                for j in range(i, len(dna)-2, 3):
                    stop = dna[j:j+3]

                    if stop in stop_codons:
                        orf_seq = dna[i:j+3]
                        orfs.append({
                            "start": i,
                            "end": j+3,
                            "sequence": orf_seq,
                            "length": len(orf_seq)
                        })
                        break
                i = j
            else:
                i += 3

    return orfs

# ================= UI =================
st.title("🧬 DNA Analyzer Pro")

dna_input = st.text_area("Enter DNA Sequence")

# ================= ANALYZE =================
if st.button("🔬 Analyze DNA"):
    dna = dna_input.upper().replace("\n", "").strip()

    if dna:
        st.session_state.dna = dna

        seq = Seq(dna)
        st.session_state.rna = str(seq.transcribe())
        st.session_state.protein = str(seq.translate())

        st.session_state.orfs = find_orfs(dna)

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

        gc = (dna.count("G") + dna.count("C")) / len(dna) * 100
        st.write(f"GC Content: {gc:.2f}%")

        # Graph
        counts = {
            "A": dna.count("A"),
            "T": dna.count("T"),
            "G": dna.count("G"),
            "C": dna.count("C")
        }

        fig, ax = plt.subplots()
        ax.bar(counts.keys(), counts.values())
        st.pyplot(fig)

        st.write(f"RNA: {st.session_state.rna}")
        st.write(f"Protein: {st.session_state.protein}")

        # ORF
        st.subheader("🧬 ORF Detection")

        if st.session_state.orfs:
            for i, orf in enumerate(st.session_state.orfs[:5]):
                st.write(f"ORF {i+1}: Start={orf['start']} End={orf['end']} Length={orf['length']}")
                st.code(orf["sequence"][:80])
        else:
            st.warning("No ORFs found")

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

        if st.session_state.blast:
            for align in st.session_state.blast.alignments[:3]:
                for hsp in align.hsps[:1]:

                    identity = (hsp.identities / hsp.align_length) * 100

                    st.markdown(f"### 🧬 {align.title}")
                    st.write(f"E-value: {hsp.expect}")
                    st.write(f"Score: {hsp.score}")
                    st.write(f"Identity: {identity:.2f}%")

# ================= NCBI =================
with tab3:
    gene = st.text_input("Search gene (e.g. BRCA1)")

    if st.button("Search NCBI"):
        if gene:
            with st.spinner("Searching NCBI..."):
                handle = Entrez.esearch(db="gene", term=gene, retmax=5)
                record = Entrez.read(handle)
                st.session_state.ncbi_results = record["IdList"]

            st.success("Results loaded!")

    if st.session_state.ncbi_results:
        for gid in st.session_state.ncbi_results:
            st.write(f"Gene ID: {gid}")

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

    # ORF
    if st.session_state.orfs:
        story.append(Spacer(1, 10))
        story.append(Paragraph("ORF Results", styles["Heading2"]))

        for i, orf in enumerate(st.session_state.orfs[:5]):
            story.append(Paragraph(f"ORF {i+1}", styles["Normal"]))
            story.append(Paragraph(f"Start: {orf['start']}", styles["Normal"]))
            story.append(Paragraph(f"End: {orf['end']}", styles["Normal"]))
            story.append(Paragraph(f"Length: {orf['length']}", styles["Normal"]))
            story.append(Spacer(1, 5))

    # NCBI
    if st.session_state.ncbi_results:
        story.append(Spacer(1, 10))
        story.append(Paragraph("NCBI Results", styles["Heading2"]))

        for gid in st.session_state.ncbi_results:
            story.append(Paragraph(f"Gene ID: {gid}", styles["Normal"]))

    # BLAST
    if st.session_state.blast:
        story.append(Spacer(1, 10))
        story.append(Paragraph("BLAST Results", styles["Heading2"]))

        for align in st.session_state.blast.alignments[:3]:
            for hsp in align.hsps[:1]:
                identity = (hsp.identities / hsp.align_length) * 100

                story.append(Paragraph(align.title, styles["Normal"]))
                story.append(Paragraph(f"E-value: {hsp.expect}", styles["Normal"]))
                story.append(Paragraph(f"Identity: {identity:.2f}%", styles["Normal"]))
                story.append(Spacer(1, 5))

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