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
                    st.success("✅ Real BLAST completed!")

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

# ------------------ PDF ------------------

def add_watermark(c):
    c.saveState()

    from reportlab.lib.pagesizes import letter
    width, height = letter

    c.setFillGray(0.92)
    c.setFont("Helvetica-Bold", 80)

    # Perfect centered diagonal watermark
    c.translate(width / 2, height / 2)
    c.rotate(45)

    c.drawCentredString(0, 0, "Smitangshu Bio Lab")

    c.restoreState()


def add_border(c):
    c.setLineWidth(2)
    c.rect(20, 20, 550, 750)

def add_page_design(canvas, doc):
    add_border(canvas)
    add_watermark(canvas)
# -------------------- PDF --------------------

def make_pdf():

    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer)

    styles = getSampleStyleSheet()
    story = []

    # -------- DNA --------
    if st.session_state.get("dna"):
        story.append(Paragraph("DNA Analysis", styles["Heading2"]))
        story.append(Paragraph(f"Length: {len(st.session_state.dna)}", styles["Normal"]))

    story.append(Spacer(1, 15))

    # -------- NCBI --------
    if st.session_state.get("gene_results"):
        story.append(Paragraph("NCBI Gene Results", styles["Heading2"]))
        story.append(Paragraph(st.session_state.gene_results[:1500], styles["Normal"]))

    story.append(Spacer(1, 15))

    # -------- BLAST --------
    if st.session_state.get("blast"):

        story.append(Paragraph("BLAST Results", styles["Heading2"]))

        for align in st.session_state.blast.alignments[:5]:
            for hsp in align.hsps[:1]:

                identity_percent = (hsp.identities / hsp.align_length) * 100

                story.append(Paragraph(f"<b>{align.title}</b>", styles["Normal"]))
                story.append(Paragraph(f"E-value: {hsp.expect}", styles["Normal"]))
                story.append(Paragraph(f"Score: {hsp.score}", styles["Normal"]))
                story.append(Paragraph(f"Identity: {identity_percent:.2f}%", styles["Normal"]))
                story.append(Paragraph(f"Alignment Length: {hsp.align_length}", styles["Normal"]))
                story.append(Spacer(1, 10))

    # -------- BUILD --------
    doc.build(story)

    buffer.seek(0)
    return buffer

# ------------------ DOWNLOAD ------------------
if st.button("Download Report"):
    pdf = make_pdf()
    st.download_button("Download PDF", pdf, "dna_report.pdf")

# ------------------ FOOTER ------------------
st.markdown("---")
st.markdown("Made with 💕 by Smitangshu")