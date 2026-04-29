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



# ---------------- TABS UI ----------------
tab1, tab2, tab3, tab4 = st.tabs([
    "🧬 DNA Analysis",
    "🧪 BLAST",
    "🔬 Gene Search",
    "📄 Report"
])

# ---------------- DNA TAB ----------------
with tab1:
    if st.session_state.get("dna"):
        dna = st.session_state.dna

        st.subheader("🧬 DNA Analysis")

        st.write(f"Length: {len(dna)}")

        if len(dna) > 0:
            gc = (dna.count("G") + dna.count("C")) / len(dna) * 100
            st.write(f"GC Content: {gc:.2f}%")

# ---------------- BLAST TAB ----------------
with tab2:
    st.subheader("🧪 BLAST")

    if st.button("Run BLAST Search"):
        if st.session_state.get("dna"):
            with st.spinner("Running BLAST..."):
                try:
                    result = NCBIWWW.qblast("blastn", "nt", st.session_state.dna)
                    blast_records = list(NCBIXML.parse(result))

                    if blast_records:
                        st.session_state.blast = blast_records[0]
                        st.success("✅ BLAST Search Completed")

                    else:
                        raise Exception("Empty BLAST")

                except Exception:
                    st.warning("⚠️ BLAST failed. Showing demo results.")

                    # -------- FAKE BLAST --------
                    class FakeHSP:
                        def __init__(self):
                            self.expect = 0.0001
                            self.score = 100
                            self.identities = 90
                            self.align_length = 100

                    class FakeAlign:
                        def __init__(self):
                            self.title = "Demo Gene Match - Homo sapiens"
                            self.hsps = [FakeHSP()]

                    class FakeBlast:
                        def __init__(self):
                            self.alignments = [FakeAlign(), FakeAlign(), FakeAlign()]

                    st.session_state.blast = FakeBlast()
        else:
            st.warning("Enter DNA first!")

    # -------- SHOW BLAST --------
    if st.session_state.get("blast"):
        st.subheader("🔍 BLAST Results")

        for align in st.session_state.blast.alignments[:3]:
            for hsp in align.hsps[:1]:

                identity = (hsp.identities / hsp.align_length) * 100

                with st.container():
                    st.markdown(f"### 🧬 {align.title}")

                    col1, col2 = st.columns(2)

                    col1.write(f"E-value: {hsp.expect}")
                    col1.write(f"Score: {getattr(hsp, 'score', 'N/A')}")

                    col2.write(f"Identity: {identity:.2f}%")
                    col2.write(f"Length: {hsp.align_length}")

                    st.divider()

# ---------------- GENE TAB ----------------
with tab3:
    st.subheader("🔬 Gene Search")

    if st.session_state.get("gene_results"):
        st.write(st.session_state.gene_results)

# ---------------- REPORT TAB ----------------
with tab4:
    st.subheader("📄 Download Report")

    if st.session_state.get("dna"):
        pdf = make_pdf()

        st.download_button(
            label="📄 Download Full Report",
            data=pdf,
            file_name="dna_report.pdf",
            mime="application/pdf"
        )

# ----------- Full PDF Section -----------
from io import BytesIO
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter


# ----------- PDF FUNCTION -----------
def make_pdf():

    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)

    styles = getSampleStyleSheet()
    story = []

    # ----------- WATERMARK -----------
    def add_watermark(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 40)
        canvas.setFillGray(0.9)
        canvas.rotate(45)
        canvas.drawString(150, -50, "GeneScope AI")  # 🔁 change text here
        canvas.restoreState()

    # ----------- BORDER -----------
    def add_border(canvas):
        canvas.setLineWidth(2)
        canvas.rect(20, 20, 550, 750)

    # ----------- DESIGN COMBINED -----------
    def add_design(canvas, doc):
        add_border(canvas)
        add_watermark(canvas, doc)

    # ================= DNA SECTION =================
    if st.session_state.get("dna"):
        dna = st.session_state.dna

        if dna and len(dna) > 0:

            # GC content
            gc = (dna.count("G") + dna.count("C")) / len(dna) * 100

            story.append(Paragraph("DNA Analysis", styles["Heading2"]))
            story.append(Spacer(1, 10))

            story.append(Paragraph(f"Sequence: {dna}", styles["Normal"]))
            story.append(Paragraph(f"Length: {len(dna)}", styles["Normal"]))
            story.append(Paragraph(f"GC Content: {gc:.2f}%", styles["Normal"]))

            # Nucleotide distribution
            counts = {
                "A": dna.count("A"),
                "T": dna.count("T"),
                "G": dna.count("G"),
                "C": dna.count("C"),
            }

            story.append(Spacer(1, 10))
            story.append(Paragraph("Nucleotide Distribution:", styles["Heading3"]))

            for base, count in counts.items():
                story.append(Paragraph(f"{base}: {count}", styles["Normal"]))

    # ================= PROTEIN =================
    if st.session_state.get("protein"):
        story.append(Spacer(1, 15))
        story.append(Paragraph("Protein Translation", styles["Heading2"]))
        story.append(Paragraph(st.session_state.protein, styles["Normal"]))

    # ================= NCBI GENE SEARCH =================
    if st.session_state.get("gene_results"):
        story.append(Spacer(1, 15))
        story.append(Paragraph("NCBI Gene Search Results", styles["Heading2"]))

        gene_text = str(st.session_state.gene_results)
        gene_text = gene_text[:1500]  # limit size for PDF

        story.append(Paragraph(gene_text, styles["Normal"]))

    # ================= BLAST =================
    if st.session_state.get("blast"):
        story.append(Spacer(1, 15))
        story.append(Paragraph("BLAST Results", styles["Heading2"]))

        try:
            for align in st.session_state.blast.alignments[:5]:
                for hsp in align.hsps[:1]:

                    identity_percent = (hsp.identities / hsp.align_length) * 100

                    story.append(Spacer(1, 10))
                    story.append(Paragraph(f"<b>{align.title}</b>", styles["Normal"]))
                    story.append(Paragraph(f"E-value: {hsp.expect}", styles["Normal"]))
                    story.append(Paragraph(f"Score: {hsp.score}", styles["Normal"]))
                    story.append(Paragraph(f"Identity: {identity_percent:.2f}%", styles["Normal"]))
                    story.append(Paragraph(f"Alignment Length: {hsp.align_length}", styles["Normal"]))

        except Exception:
            story.append(Paragraph("BLAST data not available.", styles["Normal"]))

    # ================= BUILD =================
    doc.build(
        story,
        onFirstPage=add_design,
        onLaterPages=add_design
    )

    buffer.seek(0)
    return buffer

# ----------- DOWNLOAD BUTTON -----------
if st.session_state.get("dna"):
    pdf = make_pdf()

    st.download_button(
        label="📄 Download Full Report",
        data=pdf,
        file_name="dna_report.pdf",
        mime="application/pdf"
    )

# ------------------ FOOTER ------------------
st.markdown("---")
st.markdown("Made with 💕 by Smitangshu")