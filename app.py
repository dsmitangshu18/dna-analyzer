import streamlit as st
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
from io import BytesIO
import matplotlib.pyplot as plt

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet

# ---------------- CONFIG ----------------
Entrez.email = "smitangshudas23@gmail.com"
st.set_page_config(page_title="DNA Analyzer Pro", layout="wide")

# ---------------- SESSION ----------------
for key in ["dna", "blast", "gene_results"]:
    if key not in st.session_state:
        st.session_state[key] = None

# ---------------- FUNCTIONS ----------------

def create_gc_chart(dna):
    counts = {
        "A": dna.count("A"),
        "T": dna.count("T"),
        "G": dna.count("G"),
        "C": dna.count("C"),
    }

    plt.figure(figsize=(4,3))
    plt.bar(counts.keys(), counts.values())
    plt.title("Nucleotide Distribution")
    plt.tight_layout()
    plt.savefig("chart.png")
    plt.close()


def find_orfs(dna):
    seq = Seq(dna)
    return seq.translate()


def find_motif(dna, motif="ATG"):
    positions = []
    for i in range(len(dna)):
        if dna[i:i+len(motif)] == motif:
            positions.append(i)
    return positions


def ncbi_search(query):
    handle = Entrez.esearch(db="gene", term=query, retmax=3)
    record = Entrez.read(handle)
    return record["IdList"]


def make_pdf():
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer)
    styles = getSampleStyleSheet()
    story = []

    dna = st.session_state.get("dna")

    if dna:
        story.append(Paragraph("DNA Analysis", styles["Heading2"]))
        story.append(Paragraph(f"Sequence: {dna}", styles["Normal"]))
        story.append(Paragraph(f"Length: {len(dna)}", styles["Normal"]))

        gc = (dna.count("G")+dna.count("C"))/len(dna)*100
        story.append(Paragraph(f"GC Content: {gc:.2f}%", styles["Normal"]))

        create_gc_chart(dna)
        story.append(Image("chart.png", width=300, height=200))

        protein = find_orfs(dna)
        story.append(Paragraph(f"Protein: {protein}", styles["Normal"]))

        motifs = find_motif(dna)
        story.append(Paragraph(f"Motif positions: {motifs}", styles["Normal"]))

    if st.session_state.get("gene_results"):
        story.append(Paragraph("NCBI Results", styles["Heading2"]))
        story.append(Paragraph(str(st.session_state.gene_results), styles["Normal"]))

    if st.session_state.get("blast"):
        story.append(Paragraph("BLAST Results", styles["Heading2"]))

        for align in st.session_state.blast.alignments[:3]:
            for hsp in align.hsps[:1]:
                identity = (hsp.identities / hsp.align_length) * 100

                story.append(Paragraph(align.title, styles["Normal"]))
                story.append(Paragraph(f"E-value: {hsp.expect}", styles["Normal"]))
                story.append(Paragraph(f"Identity: {identity:.2f}%", styles["Normal"]))
                story.append(Spacer(1,10))

    doc.build(story)
    buffer.seek(0)
    return buffer


# ---------------- INPUT ----------------

st.title("🧬 DNA Analyzer Pro")

dna_input = st.text_area("Enter DNA Sequence")

if dna_input:
    dna = dna_input.upper().strip()
    st.session_state.dna = dna

# ---------------- UI ----------------

st.markdown("---")

tab1, tab2, tab3, tab4 = st.tabs([
    "🧬 DNA Analysis",
    "🧪 BLAST",
    "🔬 NCBI",
    "📄 Report"
])

# ---------------- TAB 1 ----------------
with tab1:
    if st.session_state.get("dna"):
        dna = st.session_state.dna

        st.subheader("DNA Analysis")
        st.write("Length:", len(dna))

        gc = (dna.count("G")+dna.count("C"))/len(dna)*100
        st.write("GC Content:", f"{gc:.2f}%")

        create_gc_chart(dna)
        st.image("chart.png", width=300)

        st.write("Protein:", find_orfs(dna))
        st.write("Motif positions:", find_motif(dna))

# ---------------- TAB 2 ----------------
with tab2:
    if st.button("Run BLAST"):
        if st.session_state.get("dna"):
            try:
                result = NCBIWWW.qblast("blastn", "nt", st.session_state.dna)
                blast_records = list(NCBIXML.parse(result))
                st.session_state.blast = blast_records[0]
                st.success("BLAST completed")

            except:
                st.warning("BLAST failed → using fake data")

                class FakeHSP:
                    expect = 0.0001
                    identities = 90
                    align_length = 100

                class FakeAlign:
                    title = "Fake Gene Match - Homo sapiens"
                    hsps = [FakeHSP()]

                class FakeBlast:
                    alignments = [FakeAlign(), FakeAlign()]

                st.session_state.blast = FakeBlast()

    if st.session_state.get("blast"):
        for align in st.session_state.blast.alignments[:3]:
            for hsp in align.hsps[:1]:
                identity = (hsp.identities / hsp.align_length) * 100

                st.write(align.title)
                st.write("E-value:", hsp.expect)
                st.write("Identity:", f"{identity:.2f}%")

# ---------------- TAB 3 ----------------
with tab3:
    query = st.text_input("Search gene")

    if st.button("Search NCBI"):
        results = ncbi_search(query)
        st.session_state.gene_results = results

    if st.session_state.get("gene_results"):
        st.write(st.session_state.gene_results)

# ---------------- TAB 4 ----------------
with tab4:
    if st.session_state.get("dna"):
        pdf = make_pdf()

        st.download_button(
            "Download Report",
            data=pdf,
            file_name="dna_report.pdf",
            mime="application/pdf"
        )

# ------------------ FOOTER ------------------
st.markdown("---")
st.markdown("Made with 💕 by Smitangshu")