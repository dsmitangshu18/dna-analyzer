import streamlit as st
from Bio.Seq import Seq
from Bio import Entrez
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.pdfgen import canvas
from reportlab.lib import colors
import io

# ---------------- CONFIG ----------------
Entrez.email = "smitangshudas23@gmail.com"

st.set_page_config(page_title="DNA Analyzer", layout="wide")

st.title("🧬 DNA Sequence Analyzer")
st.caption("Professional Bioinformatics Tool")

# ---------------- INPUT ----------------
dna_input = st.text_area("Enter DNA Sequence")

# ---------------- SESSION ----------------
if "analysis_done" not in st.session_state:
    st.session_state.analysis_done = False

if "ncbi_results" not in st.session_state:
    st.session_state.ncbi_results = ""

if "analysis_data" not in st.session_state:
    st.session_state.analysis_data = {}

# ---------------- ANALYZE ----------------
if st.button("Analyze"):
    if dna_input:
        dna = dna_input.upper()
        if set(dna).issubset({"A","T","G","C"}):
            st.session_state.analysis_done = True
            st.session_state.dna = dna

            seq = Seq(dna)

            reverse_comp = str(seq.reverse_complement())
            rna = str(seq.transcribe())
            protein = str(seq.translate())

            length = len(dna)

            gc = (dna.count("G") + dna.count("C")) / length * 100

            counts = {
                "A": dna.count("A"),
                "T": dna.count("T"),
                "G": dna.count("G"),
                "C": dna.count("C")
            }

            percentages = {
                "A": (counts["A"] / length) * 100,
                "T": (counts["T"] / length) * 100,
                "G": (counts["G"] / length) * 100,
                "C": (counts["C"] / length) * 100,
            }

            st.session_state.analysis_data = {
                "dna": dna,
                "reverse": reverse_comp,
                "rna": rna,
                "protein": protein,
                "gc": gc,
                "length": length,
                "counts": counts,
                "percentages": percentages
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
        f"A: {data['percentages']['A']:.2f}% | "
        f"T: {data['percentages']['T']:.2f}% | "
        f"G: {data['percentages']['G']:.2f}% | "
        f"C: {data['percentages']['C']:.2f}%"
    )

    fig, ax = plt.subplots()
    ax.bar(data["counts"].keys(), data["counts"].values())
    st.pyplot(fig)

# ---------------- NCBI SEARCH ----------------
st.subheader("🌐 NCBI Sequence Search")

ncbi_query = st.text_input("Enter DNA / Gene for NCBI search")

if st.button("Search NCBI Database"):
    if not ncbi_query:
        st.warning("Enter something to search")
    else:
        try:
            with st.spinner("Searching NCBI..."):

                handle = Entrez.esearch(
                    db="nucleotide",
                    term=ncbi_query,
                    retmax=3
                )
                record = Entrez.read(handle)
                handle.close()

                ids = record["IdList"]

                if not ids:
                    st.warning("No matches found")
                else:
                    st.success(f"Found {len(ids)} matches")

                    fetch = Entrez.efetch(
                        db="nucleotide",
                        id=",".join(ids),
                        rettype="fasta",
                        retmode="text"
                    )

                    results = fetch.read()
                    fetch.close()

                    st.session_state.ncbi_results = results

        except Exception as e:
            st.error(f"Error: {e}")

# ---------------- SHOW NCBI ----------------
if st.session_state.ncbi_results:
    st.subheader("📄 Retrieved sequences")
    st.text_area("Preview", st.session_state.ncbi_results[:2000], height=300)

# ---------------- PDF WATERMARK + BORDER ----------------
def add_watermark_and_border(c, doc):
    width, height = A4

    # watermark
    c.setFont("Helvetica", 40)
    c.setFillColorRGB(0.9, 0.9, 0.9)
    c.saveState()
    c.translate(width/2, height/2)
    c.rotate(45)
    c.drawCentredString(0, 0, "SMITANGSHU BIO LAB")
    c.restoreState()

    # border
    c.setStrokeColor(colors.black)
    c.rect(30, 30, width-60, height-60)

# ---------------- CREATE PDF ----------------
def create_pdf():
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4)
    styles = getSampleStyleSheet()

    content = []

    content.append(Paragraph("DNA Analysis Report", styles["Title"]))
    content.append(Spacer(1, 10))

    if st.session_state.analysis_done:
        data = st.session_state.analysis_data

        content.append(Paragraph(f"<b>DNA:</b> {data['dna']}", styles["Normal"]))
        content.append(Paragraph(f"<b>Reverse:</b> {data['reverse']}", styles["Normal"]))
        content.append(Paragraph(f"<b>RNA:</b> {data['rna']}", styles["Normal"]))
        content.append(Paragraph(f"<b>Protein:</b> {data['protein']}", styles["Normal"]))

        content.append(Paragraph(f"<b>Length:</b> {data['length']} bases", styles["Normal"]))
        content.append(Paragraph(f"<b>GC Content:</b> {data['gc']:.2f}%", styles["Normal"]))

        counts = data["counts"]
        content.append(Paragraph(
            f"<b>Counts:</b> A:{counts['A']} T:{counts['T']} G:{counts['G']} C:{counts['C']}",
            styles["Normal"]
        ))

        perc = data["percentages"]
        content.append(Paragraph(
            f"<b>Percentages:</b> A:{perc['A']:.2f}% T:{perc['T']:.2f}% G:{perc['G']:.2f}% C:{perc['C']:.2f}%",
            styles["Normal"]
        ))

        content.append(Spacer(1, 10))

    if st.session_state.ncbi_results:
        content.append(Paragraph("<b>NCBI Results:</b>", styles["Heading2"]))
        content.append(Paragraph(st.session_state.ncbi_results[:1000], styles["Normal"]))

    doc.build(content, onFirstPage=add_watermark_and_border, onLaterPages=add_watermark_and_border)

    buffer.seek(0)
    return buffer

# ---------------- DOWNLOAD ----------------
if st.button("Download Professional Report"):
    pdf = create_pdf()
    st.download_button(
        label="Download PDF",
        data=pdf,
        file_name="dna_report.pdf",
        mime="application/pdf"
    )

# ---------------- FOOTER ----------------
st.markdown("---")
st.markdown("❤️ Made with 💕 by Smitangshu")