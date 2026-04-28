import streamlit as st
from Bio.Seq import Seq
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from io import BytesIO

Entrez.email = "your_email@example.com"

st.set_page_config(page_title="DNA Analyzer", layout="wide")

# ---------------- STATE ----------------
if "analysis" not in st.session_state:
    st.session_state.analysis = None
if "ncbi" not in st.session_state:
    st.session_state.ncbi = None
if "blast" not in st.session_state:
    st.session_state.blast = None

# ---------------- HEADER ----------------
st.title("🧬 DNA Analyzer")

dna_input = st.text_area("Enter DNA Sequence")

def valid_dna(seq):
    return all(c in "ATGC" for c in seq.upper())

# ---------------- ANALYZE ----------------
if st.button("Analyze"):
    if not dna_input or not valid_dna(dna_input):
        st.error("Invalid DNA sequence")
    else:
        dna = dna_input.upper()
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

    st.write("Nucleotide Counts:", data["counts"])
    st.write("RNA:", data["rna"])
    st.write("Protein:", data["protein"])

    # Graph
    fig, ax = plt.subplots()
    ax.bar(data["counts"].keys(), data["counts"].values())
    st.pyplot(fig)

# ---------------- NCBI SEARCH ----------------
st.subheader("🌐 NCBI Sequence Search")

if st.button("Search NCBI"):
    if st.session_state.analysis:
        dna = st.session_state.analysis["dna"]

        with st.spinner("Searching NCBI..."):
            handle = Entrez.esearch(
                db="nucleotide",
                term=f"{dna[:50]}[All Fields]",
                retmax=3
            )
            record = Entrez.read(handle)
            ids = record["IdList"]

            if ids:
                fetch = Entrez.efetch(
                    db="nucleotide",
                    id=",".join(ids),
                    rettype="fasta",
                    retmode="text"
                )
                results = fetch.read()
                st.session_state.ncbi = results
            else:
                st.warning("No matches found")

# SHOW NCBI
if st.session_state.ncbi:
    st.subheader("📄 NCBI Results")
    st.text_area("Preview", st.session_state.ncbi[:2000], height=300)

# ---------------- GENE SEARCH ----------------
st.subheader("🧬 NCBI Gene Search")
gene = st.text_input("Enter gene name")

if st.button("Search Gene"):
    if gene:
        handle = Entrez.esearch(db="nucleotide", term=gene, retmax=2)
        record = Entrez.read(handle)
        ids = record["IdList"]

        if ids:
            fetch = Entrez.efetch(db="nucleotide", id=",".join(ids),
                                 rettype="fasta", retmode="text")
            results = fetch.read()
            st.text_area("Gene Results", results[:2000], height=300)

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
        dna = st.session_state.analysis["dna"]

        with st.spinner("Running BLAST..."):
            result = NCBIWWW.qblast("blastn", "nt", dna)
            blast_record = NCBIXML.read(result)

            st.session_state.blast = blast_record

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

            with st.expander(f"{align.title[:60]}"):
                st.write(f"Score: {hsp.score}")
                st.write(f"E-value: {hsp.expect}")
                st.write(f"Quality: {quality(hsp.expect)}")

                st.code(f"""
Query: {hsp.query}
Match: {hsp.match}
Sbjct: {hsp.sbjct}
""")

# ---------------- PDF ----------------
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
        story.append(Paragraph(f"Counts: {a['counts']}", styles['Normal']))
        story.append(Paragraph(f"RNA: {a['rna']}", styles['Normal']))
        story.append(Paragraph(f"Protein: {a['protein']}", styles['Normal']))

    if st.session_state.ncbi:
        story.append(Paragraph("NCBI Results", styles['Heading2']))
        story.append(Paragraph(st.session_state.ncbi[:1000], styles['Normal']))

    if st.session_state.blast:
        story.append(Paragraph("BLAST Results", styles['Heading2']))
        for align in st.session_state.blast.alignments[:2]:
            for hsp in align.hsps:
                story.append(Paragraph(f"{align.title}", styles['Normal']))
                story.append(Paragraph(f"E-value: {hsp.expect}", styles['Normal']))

    story.append(Spacer(1, 20))
    story.append(Paragraph("Made with 💕 by Smitangshu", styles['Normal']))

    doc.build(story)
    buffer.seek(0)
    return buffer

if st.button("Download Report"):
    pdf = make_pdf()
    st.download_button("Download PDF", pdf, "report.pdf")

# ---------------- FOOTER ----------------
st.markdown("---")
st.markdown("Made with 💕 by Smitangshu")