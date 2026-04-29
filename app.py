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

# -------------------- BLAST --------------------
st.subheader("🔬 BLAST")

if st.button("Run BLAST Search"):
    if st.session_state.get("dna"):

        with st.spinner("Running BLAST..."):
            try:
                from Bio.Blast import NCBIWWW, NCBIXML

                result = NCBIWWW.qblast("blastn", "nt", st.session_state.dna)
                blast_records = list(NCBIXML.parse(result))

                if blast_records and len(blast_records) > 0:
                    st.session_state.blast = blast_records[0]
                    st.success("🧬Sequence alignment completed!")

                else:
                    raise Exception("Empty BLAST result")

            except Exception as e:
                st.warning("⚠️ Real BLAST failed. Showing demo results instead.")

                # -------- FAKE BLAST (FALLBACK) --------
                class FakeHSP:
                    def __init__(self):
                        self.expect = 0.0001
                        self.score = 120
                        self.align_length = 85
                        self.identities = 75
                        self.query = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
                        self.match = "||||||||||||||||||||||||||||||||||||||||"
                        self.sbjct = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"

                class FakeAlign:
                    def __init__(self):
                        self.title = "Fake Gene Match - Homo sapiens"
                        self.hsps = [FakeHSP()]

                class FakeBlast:
                    def __init__(self):
                        self.alignments = [FakeAlign(), FakeAlign(), FakeAlign()]

                st.session_state.blast = FakeBlast()

    else:
        st.warning("⚠️ Run DNA analysis first!")

# -------------------- SHOW BLAST --------------------
if st.session_state.get("blast"):

    st.subheader("🔎 BLAST Results")

    best_match = None
    best_evalue = float("inf")

    found = False

    for align in st.session_state.blast.alignments[:5]:
        for hsp in align.hsps[:1]:

            found = True

            # Calculate identity %
            identity_percent = (hsp.identities / hsp.align_length) * 100

            # Track BEST match (lowest E-value)
            if hsp.expect < best_evalue:
                best_evalue = hsp.expect
                best_match = align.title

            # -------- UI BLOCK --------
            with st.container():
                st.markdown(f"### 🧬 {align.title}")

                col1, col2 = st.columns(2)

                col1.write(f"**E-value:** {hsp.expect}")
                col1.write(f"**Score:** {getattr(hsp, 'score', 'N/A')}")

                col2.write(f"**Identity:** {identity_percent:.2f}%")
                col2.write(f"**Length:** {hsp.align_length}")

                st.divider()

    # If nothing found
    if not found:
        st.warning("No alignments found.")

    # Show BEST MATCH
    if best_match:
        st.success(f"🏆 Best Match: {best_match} (E-value: {best_evalue})")

# -------------------- PDF SECTION --------------------

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from io import BytesIO

# ---------- DESIGN FUNCTIONS ----------

def add_border(c):
    c.setLineWidth(2)
    c.rect(20, 20, 550, 750)

def add_watermark(c):
    c.saveState()
    c.setFont("Helvetica", 40)
    c.setFillGray(0.9)
    c.drawCentredString(300, 400, "DNA Analyzer")
    c.restoreState()

def add_page_design(c, doc):
    add_border(c)
    add_watermark(c)

# ---------- PDF GENERATOR ----------

def make_pdf():
    buffer = BytesIO()

    doc = SimpleDocTemplate(buffer, pagesize=letter)
    styles = getSampleStyleSheet()
    story = []

    # -------- DNA SECTION --------
    if st.session_state.get("dna"):

        dna = st.session_state.dna

        story.append(Paragraph("DNA Analysis", styles["Heading2"]))
        story.append(Paragraph(f"Length: {len(dna)}", styles["Normal"]))
        story.append(Spacer(1, 10))

        story.append(Paragraph("<b>Sequence:</b>", styles["Normal"]))
        story.append(Paragraph(dna[:1000], styles["Normal"]))
        story.append(Spacer(1, 10))

        if st.session_state.get("gc"):
            story.append(Paragraph(f"GC Content: {st.session_state.gc}%", styles["Normal"]))

        story.append(Spacer(1, 20))

    # -------- BLAST SECTION --------
    if st.session_state.get("blast"):

        story.append(Paragraph("BLAST Results", styles["Heading2"]))
        story.append(Spacer(1, 10))

        best_match = None
        best_evalue = float("inf")

        for align in st.session_state.blast.alignments[:5]:
            for hsp in align.hsps[:1]:

                identity_percent = (hsp.identities / hsp.align_length) * 100

                # Track best match
                if hsp.expect < best_evalue:
                    best_evalue = hsp.expect
                    best_match = align.title

                story.append(Paragraph(f"<b>{align.title}</b>", styles["Normal"]))
                story.append(Paragraph(f"E-value: {hsp.expect}", styles["Normal"]))
                story.append(Paragraph(f"Score: {getattr(hsp, 'score', 'N/A')}", styles["Normal"]))
                story.append(Paragraph(f"Identity: {identity_percent:.2f}%", styles["Normal"]))
                story.append(Paragraph(f"Alignment Length: {hsp.align_length}", styles["Normal"]))
                story.append(Spacer(1, 12))

        # -------- BEST MATCH --------
        if best_match:
            story.append(Spacer(1, 10))
            story.append(Paragraph("<b>Best Match:</b>", styles["Heading3"]))
            story.append(Paragraph(best_match, styles["Normal"]))
            story.append(Paragraph(f"E-value: {best_evalue}", styles["Normal"]))

    # -------- BUILD PDF --------
    doc.build(
        story,
        onFirstPage=add_page_design,
        onLaterPages=add_page_design
    )

    buffer.seek(0)
    return buffer

# ---------- DOWNLOAD BUTTON ----------

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