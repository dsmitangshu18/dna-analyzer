import streamlit as st
from Bio.Seq import Seq
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
import io

# ---------------- CONFIG ----------------
Entrez.email = "smitangshudas23@gmail.com"
st.set_page_config(page_title="DNA Analyzer", layout="wide")

st.title("🧬 DNA Sequence Analyzer")
st.caption("Professional Bioinformatics Tool")

# ---------------- SESSION ----------------
if "analysis_done" not in st.session_state:
    st.session_state.analysis_done = False

if "analysis_data" not in st.session_state:
    st.session_state.analysis_data = {}

if "ncbi_results" not in st.session_state:
    st.session_state.ncbi_results = ""

if "blast_results" not in st.session_state:
    st.session_state.blast_results = None

# ---------------- INPUT ----------------
dna_input = st.text_area("Enter DNA Sequence")

# ---------------- ANALYZE ----------------
if st.button("Analyze"):
    if dna_input:
        dna = dna_input.upper()
        if set(dna).issubset({"A","T","G","C"}):
            seq = Seq(dna)

            reverse = str(seq.reverse_complement())
            rna = str(seq.transcribe())
            protein = str(seq.translate())

            length = len(dna)
            gc = (dna.count("G") + dna.count("C")) / length * 100

            counts = {b: dna.count(b) for b in "ATGC"}
            perc = {b: (counts[b]/length)*100 for b in "ATGC"}

            st.session_state.analysis_done = True
            st.session_state.analysis_data = {
                "dna": dna,
                "reverse": reverse,
                "rna": rna,
                "protein": protein,
                "length": length,
                "gc": gc,
                "counts": counts,
                "perc": perc
            }
        else:
            st.error("Invalid DNA sequence")
    else:
        st.warning("Enter DNA first")

# ---------------- SHOW ANALYSIS ----------------
if st.session_state.analysis_done:
    data = st.session_state.analysis_data

    st.subheader("🧪 Sequences")
    st.code(data["dna"])
    st.code(data["reverse"])
    st.code(data["rna"])
    st.code(data["protein"])

    st.success(f"GC Content: {data['gc']:.2f}%")
    st.info(f"Length: {data['length']} bases")

    st.write("### Nucleotide Percentages")
    st.write(
        f"A: {data['perc']['A']:.2f}% | "
        f"T: {data['perc']['T']:.2f}% | "
        f"G: {data['perc']['G']:.2f}% | "
        f"C: {data['perc']['C']:.2f}%"
    )

    fig, ax = plt.subplots()
    ax.bar(data["counts"].keys(), data["counts"].values())
    st.pyplot(fig)

    # -------- MOTIF --------
    st.subheader("🔍 Motif Finder")
    motifs = st.text_input("Motifs (comma separated)", "ATG")
    for m in [x.strip().upper() for x in motifs.split(",") if x]:
        pos = [i for i in range(len(data["dna"])) if data["dna"].startswith(m, i)]
        st.success(f"{m} → {pos}")

    # -------- ORF --------
    st.subheader("🧬 ORF Finder")
    dna = data["dna"]
    stops = ["TAA","TAG","TGA"]
    for i in range(len(dna)-2):
        if dna[i:i+3]=="ATG":
            for j in range(i, len(dna)-2, 3):
                if dna[j:j+3] in stops:
                    st.code(dna[i:j+3])
                    break

# ---------------- NCBI ----------------
st.subheader("🌐 NCBI Search")
query = st.text_input("Enter DNA / Gene")

if st.button("Search NCBI"):
    if query:
        try:
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=2)
            ids = Entrez.read(handle)["IdList"]
            handle.close()

            if ids:
                fetch = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text")
                st.session_state.ncbi_results = fetch.read()
                fetch.close()
            else:
                st.warning("No matches")
        except Exception as e:
            st.error(e)

if st.session_state.ncbi_results:
    st.subheader("📄 NCBI Results")
    st.text_area("Preview", st.session_state.ncbi_results[:2000], height=300)

# ---------------- BLAST ----------------
st.subheader("🧬 BLAST Search")

if st.button("Run BLAST"):
    if st.session_state.analysis_done:
        dna = st.session_state.analysis_data["dna"]
        try:
            result = NCBIWWW.qblast("blastn","nt",dna)
            record = NCBIXML.read(result)
            st.session_state.blast_results = record
        except Exception as e:
            st.error(e)
    else:
        st.warning("Analyze first")

if st.session_state.blast_results:
    st.subheader("🔬 BLAST Results")
    for align in st.session_state.blast_results.alignments[:2]:
        st.markdown(f"**{align.title}**")
        for hsp in align.hsps[:1]:
            st.write(f"Score: {hsp.score} | E-value: {hsp.expect}")
            st.code(hsp.query[:100])
            st.code(hsp.match[:100])
            st.code(hsp.sbjct[:100])

# ---------------- PDF ----------------
from reportlab.pdfgen import canvas

def decorate(c, doc):
    w,h = A4
    c.setFont("Helvetica",40)
    c.setFillColorRGB(0.9,0.9,0.9)
    c.saveState()
    c.translate(w/2,h/2)
    c.rotate(45)
    c.drawCentredString(0,0,"DNA ANALYZER")
    c.restoreState()
    c.rect(30,30,w-60,h-60)

def create_pdf():
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4)
    styles = getSampleStyleSheet()
    content = []

    content.append(Paragraph("DNA Report", styles["Title"]))
    content.append(Spacer(1,10))

    if st.session_state.analysis_done:
        d = st.session_state.analysis_data
        content.append(Paragraph(f"DNA: {d['dna']}", styles["Normal"]))
        content.append(Paragraph(f"GC: {d['gc']:.2f}%", styles["Normal"]))
        content.append(Paragraph(f"Length: {d['length']}", styles["Normal"]))

    if st.session_state.ncbi_results:
        content.append(Paragraph("NCBI Results:", styles["Heading2"]))
        content.append(Paragraph(st.session_state.ncbi_results[:800], styles["Normal"]))

    if st.session_state.blast_results:
        content.append(Paragraph("BLAST Results:", styles["Heading2"]))
        for a in st.session_state.blast_results.alignments[:1]:
            content.append(Paragraph(a.title, styles["Normal"]))

    doc.build(content, onFirstPage=decorate, onLaterPages=decorate)
    buffer.seek(0)
    return buffer

if st.button("Download Report"):
    pdf = create_pdf()
    st.download_button("Download PDF", pdf, "dna_report.pdf")

# ---------------- FOOTER ----------------
st.markdown("---")
st.markdown("❤️ Made with 💕 by Smitangshu")