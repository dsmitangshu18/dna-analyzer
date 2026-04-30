import streamlit as st
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import requests
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from io import BytesIO

# ---------------- PAGE CONFIG ----------------
st.set_page_config(page_title="DNA Analyzer Pro", layout="wide")

# ---------------- HELPERS ----------------

def validate_dna(seq):
    return all(c in "ATGC" for c in seq.upper())

def get_gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) * 100

def nucleotide_counts(seq):
    return {
        "A": seq.count("A"),
        "T": seq.count("T"),
        "G": seq.count("G"),
        "C": seq.count("C")
    }

def translate_dna(seq):
    return str(Seq(seq).translate())

def transcribe_dna(seq):
    return str(Seq(seq).transcribe())

def find_orfs(seq):
    orfs = []
    for i in range(len(seq)-2):
        if seq[i:i+3] == "ATG":
            for j in range(i, len(seq)-2, 3):
                codon = seq[j:j+3]
                if codon in ["TAA","TAG","TGA"]:
                    orfs.append((i, j+3))
                    break
    return orfs

def fake_blast(seq):
    return [
        {"title": "Homo sapiens hemoglobin subunit beta",
         "score": 200,
         "identity": "99%",
         "evalue": "1e-50"}
    ]

def ncbi_search(gene):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene}&retmode=json"
    try:
        res = requests.get(url).json()
        return res["esearchresult"]["idlist"]
    except:
        return []

# ---------------- PDF ----------------

def add_design(canvas, doc):
    canvas.setLineWidth(2)
    canvas.rect(20, 20, 550, 750)
    canvas.setFont("Helvetica", 40)
    canvas.setFillGray(0.9)
    canvas.drawCentredString(300, 400, "DNA ANALYZER PRO")

def make_pdf(data):
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)
    styles = getSampleStyleSheet()
    story = []

    dna = data["dna"]

    story.append(Paragraph("DNA Analysis", styles["Heading2"]))
    story.append(Spacer(1,10))
    story.append(Paragraph(f"Sequence: {dna}", styles["Normal"]))
    story.append(Paragraph(f"Length: {len(dna)}", styles["Normal"]))
    story.append(Paragraph(f"GC Content: {data['gc']:.2f}%", styles["Normal"]))
    story.append(Paragraph(f"RNA: {data['rna']}", styles["Normal"]))
    story.append(Paragraph(f"Protein: {data['protein']}", styles["Normal"]))

    story.append(Spacer(1,15))
    story.append(Paragraph("ORF Results", styles["Heading2"]))
    for start,end in data["orfs"]:
        story.append(Paragraph(f"Start: {start}, End: {end}", styles["Normal"]))

    story.append(Spacer(1,15))
    story.append(Paragraph("BLAST Results", styles["Heading2"]))
    for r in data["blast"]:
        story.append(Paragraph(r["title"], styles["Normal"]))

    story.append(Spacer(1,15))
    story.append(Paragraph("NCBI Results", styles["Heading2"]))
    for g in data["ncbi"]:
        story.append(Paragraph(f"Gene ID: {g}", styles["Normal"]))

    doc.build(story, onFirstPage=add_design, onLaterPages=add_design)
    buffer.seek(0)
    return buffer

# ---------------- UI ----------------

st.title("🧬 DNA Analyzer Pro")

dna_input = st.text_area("Enter DNA Sequence")

if st.button("🔍 Analyze Sequence"):

    if not dna_input:
        st.warning("Enter DNA sequence first")
    elif not validate_dna(dna_input):
        st.error("Invalid DNA sequence")
    else:
        dna = dna_input.upper()

        gc = get_gc_content(dna)
        counts = nucleotide_counts(dna)
        protein = translate_dna(dna)
        rna = transcribe_dna(dna)
        orfs = find_orfs(dna)

        st.session_state.data = {
            "dna": dna,
            "gc": gc,
            "counts": counts,
            "protein": protein,
            "rna": rna,
            "orfs": orfs,
            "blast": [],
            "ncbi": []
        }

# ---------------- TABS ----------------

tab1, tab2, tab3, tab4 = st.tabs(["DNA", "BLAST", "NCBI", "Report"])

# DNA TAB
with tab1:
    if "data" in st.session_state:
        d = st.session_state.data

        st.metric("Length", len(d["dna"]))
        st.metric("GC %", f"{d['gc']:.2f}")

        fig, ax = plt.subplots()
        ax.bar(d["counts"].keys(), d["counts"].values())
        st.pyplot(fig)

        st.code(d["protein"], language="text")
        st.code(d["rna"], language="text")

        st.subheader("ORFs")
        for o in d["orfs"]:
            st.write(o)

# BLAST TAB
with tab2:
    if "data" in st.session_state:
        if st.button("Run BLAST"):
            with st.spinner("Running BLAST..."):
                st.session_state.data["blast"] = fake_blast(st.session_state.data["dna"])

        for r in st.session_state.data["blast"]:
            st.write(r["title"])

# NCBI TAB
with tab3:
    gene = st.text_input("Search gene")

    if st.button("Search NCBI"):
        with st.spinner("Searching..."):
            results = ncbi_search(gene)
            st.session_state.data["ncbi"] = results

    if "data" in st.session_state:
        for g in st.session_state.data["ncbi"]:
            st.write(f"https://www.ncbi.nlm.nih.gov/gene/{g}")

# REPORT TAB
with tab4:
    if "data" in st.session_state:
        pdf = make_pdf(st.session_state.data)

        st.download_button(
            "Download PDF",
            data=pdf,
            file_name="dna_report.pdf",
            mime="application/pdf"
        )

# ------------------ FOOTER ------------------
st.markdown("---")
st.markdown("Made with 💕 by Smitangshu")