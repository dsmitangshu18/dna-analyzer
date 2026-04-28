import streamlit as st
from Bio.Seq import Seq
from Bio import Entrez
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet

# ---------------- CONFIG ----------------
st.set_page_config(page_title="DNA Analyzer", layout="wide")
Entrez.email = "smitangshudas23@gmail.com"

# ---------------- SESSION ----------------
if "analysis" not in st.session_state:
    st.session_state.analysis = None

if "ncbi_result" not in st.session_state:
    st.session_state.ncbi_result = None

# ---------------- TITLE ----------------
st.title("🧬 DNA Sequence Analyzer")

# ---------------- INPUT ----------------
dna_input = st.text_area("Enter DNA Sequence")

# ---------------- VALIDATION ----------------
def is_valid_dna(seq):
    return all(base in "ATCG" for base in seq.upper())

# ---------------- ANALYZE ----------------
if st.button("Analyze"):

    if not dna_input:
        st.warning("Enter DNA first")
    elif not is_valid_dna(dna_input):
        st.error("Invalid DNA sequence")
    else:
        dna = dna_input.upper()
        seq = Seq(dna)

        reverse_comp = str(seq.reverse_complement())
        rna = str(seq.transcribe())
        protein = str(seq.translate())

        # GC CONTENT
        gc = (dna.count("G") + dna.count("C")) / len(dna) * 100

        # NUCLEOTIDE COUNT
        counts = {
            "A": dna.count("A"),
            "T": dna.count("T"),
            "C": dna.count("C"),
            "G": dna.count("G"),
        }

        # SAVE IN SESSION
        st.session_state.analysis = {
            "dna": dna,
            "reverse": reverse_comp,
            "rna": rna,
            "protein": protein,
            "gc": gc,
            "counts": counts
        }

# ---------------- SHOW ANALYSIS ----------------
if st.session_state.analysis:

    data = st.session_state.analysis

    st.subheader("🧪 Sequences")
    st.code(data["dna"])
    st.code(data["reverse"])
    st.code(data["rna"])
    st.code(data["protein"])

    st.subheader("📊 GC Content")
    st.success(f"GC Content: {data['gc']:.2f}%")

    st.subheader("📈 Nucleotide Distribution")

    fig, ax = plt.subplots()
    ax.bar(data["counts"].keys(), data["counts"].values())
    st.pyplot(fig)

    # -------- MOTIF --------
    st.subheader("🔍 Motif Finder")
    motifs_input = st.text_input("Enter motifs", "ATG")
    motifs = [m.strip().upper() for m in motifs_input.split(",") if m.strip()]

    for motif in motifs:
        positions = [i for i in range(len(data["dna"])) if data["dna"].startswith(motif, i)]
        st.success(f"{motif} → {positions}")

# ---------------- NCBI SEARCH ----------------
st.subheader("🌐 NCBI Sequence Search")

if st.button("Search NCBI"):

    if not st.session_state.analysis:
        st.warning("Run Analyze first")
    else:
        try:
            dna = st.session_state.analysis["dna"]

            with st.spinner("Searching NCBI..."):
                handle = Entrez.esearch(
                    db="nucleotide",
                    term=dna[:50],
                    retmax=2
                )
                record = Entrez.read(handle)
                handle.close()

                ids = record["IdList"]

                if not ids:
                    st.warning("No matches found")
                else:
                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        id=",".join(ids),
                        rettype="fasta",
                        retmode="text"
                    )

                    results = fetch_handle.read()
                    fetch_handle.close()

                    st.session_state.ncbi_result = results

        except Exception as e:
            st.error(f"Error: {e}")

# ---------------- SHOW NCBI ----------------
if st.session_state.ncbi_result:

    st.subheader("📄 NCBI Results")
    st.code(st.session_state.ncbi_result[:1000])

# ---------------- PDF DOWNLOAD ----------------
st.subheader("📥 Download Full Report")

if st.button("Download PDF"):

    if not st.session_state.analysis:
        st.warning("Analyze first")
    else:
        doc = SimpleDocTemplate("report.pdf")
        styles = getSampleStyleSheet()

        elements = []

        data = st.session_state.analysis

        elements.append(Paragraph("DNA Analysis Report", styles["Title"]))
        elements.append(Spacer(1, 10))

        elements.append(Paragraph(f"DNA: {data['dna']}", styles["Normal"]))
        elements.append(Paragraph(f"GC Content: {data['gc']:.2f}%", styles["Normal"]))
        elements.append(Paragraph(f"Protein: {data['protein']}", styles["Normal"]))

        if st.session_state.ncbi_result:
            elements.append(Spacer(1, 10))
            elements.append(Paragraph("NCBI Results:", styles["Heading2"]))
            elements.append(Paragraph(st.session_state.ncbi_result[:1000], styles["Normal"]))

        doc.build(elements)

        with open("report.pdf", "rb") as f:
            st.download_button(
                "⬇️ Download Report",
                f,
                file_name="dna_report.pdf"
            )