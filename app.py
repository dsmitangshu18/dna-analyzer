import streamlit as st
from Bio.Seq import Seq
import pandas as pd
from datetime import datetime

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.pdfgen import canvas

# ---------------- PAGE ----------------
st.set_page_config(page_title="DNA Analyzer", layout="wide")

st.title("🧬 DNA Sequence Analyzer")
st.caption("Professional Bioinformatics Tool")

# ---------------- SESSION ----------------
if "dna" not in st.session_state:
    st.session_state["dna"] = None

# ---------------- INPUT ----------------
uploaded_file = st.file_uploader("Upload FASTA file", type=["txt", "fasta", "fa"])
dna_text = st.text_area("Enter DNA Sequence")

# ---------------- ANALYZE ----------------
if st.button("Analyze"):

    if uploaded_file:
        content = uploaded_file.read().decode("utf-8")
        dna = "".join([l for l in content.split("\n") if not l.startswith(">")])
    else:
        dna = dna_text

    dna = dna.upper().replace(" ", "").replace("\n", "")

    if not dna:
        st.warning("Enter DNA sequence")
        st.stop()

    if not set(dna).issubset(set("ATGC")):
        st.error("Invalid DNA sequence")
        st.stop()

    st.session_state["dna"] = dna


# ---------------- MAIN ----------------
if st.session_state["dna"]:

    dna = st.session_state["dna"]
    sequence = Seq(dna)

    gc = (dna.count("G") + dna.count("C")) / len(dna) * 100
    at = 100 - gc

    reverse_comp = str(sequence.reverse_complement())
    rna = str(sequence.transcribe())
    protein = str(sequence.translate(to_stop=True))

    # ---------------- RESULTS ----------------
    st.subheader("📊 Results")
    col1, col2 = st.columns(2)
    col1.metric("Length", len(dna))
    col2.metric("GC %", round(gc, 2))

    # ---------------- NUCLEOTIDE GRAPH ----------------
    st.subheader("🧪 Nucleotide Count")

    counts = {
        "A": dna.count("A"),
        "T": dna.count("T"),
        "G": dna.count("G"),
        "C": dna.count("C")
    }
    st.bar_chart(counts)

    # ---------------- GC GRAPH ----------------
    st.subheader("📈 GC Content Across Sequence")

    window_size = st.slider("GC Window Size", 3, 20, 5)

    gc_values = []
    for i in range(len(dna) - window_size + 1):
        window = dna[i:i+window_size]
        gc_percent = (window.count("G") + window.count("C")) / window_size * 100
        gc_values.append(gc_percent)

    st.line_chart(gc_values)

    # ---------------- SEQUENCES ----------------
    st.subheader("🧬 Sequences")
    st.code(dna)
    st.code(reverse_comp)
    st.code(rna)
    st.code(protein)

    # ---------------- MOTIF ----------------
    st.subheader("🔍 Motif Finder")
    motifs_input = st.text_input("Enter motifs (comma separated)", "ATG")
    motifs = [m.strip().upper() for m in motifs_input.split(",") if m.strip()]

    for motif in motifs:
        positions = [i for i in range(len(dna)) if dna.startswith(motif, i)]
        st.success(f"{motif} → Positions: {positions}")

    # ---------------- ORF ----------------
    st.subheader("🧬 ORF Finder")

    def find_orfs(dna):
        start = "ATG"
        stops = ["TAA", "TAG", "TGA"]
        orfs = []

        for frame in range(3):
            i = frame
            while i < len(dna) - 2:
                if dna[i:i+3] == start:
                    for j in range(i, len(dna)-2, 3):
                        if dna[j:j+3] in stops:
                            seq = dna[i:j+3]
                            orfs.append((frame, i, j+3, seq))
                            break
                i += 3
        return orfs

    orfs = find_orfs(dna)

    if orfs:
        for i, (frame, start, end, seq) in enumerate(orfs, 1):
            with st.expander(f"ORF {i} (Length {len(seq)})"):
                  st.write(f"Frame: {frame} | Start: {start} | End: {end}")
                  st.code(seq)
    else:
         st.info("No ORFs found in this sequence")

    # ---------------- PDF ----------------
    def draw_border(canvas, doc):
        canvas.setLineWidth(2)
        canvas.rect(20, 20, 555, 802)

    def draw_watermark(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica-Bold", 50)
        canvas.setFillColorRGB(0.9, 0.9, 0.9)
        canvas.translate(300, 400)
        canvas.rotate(45)
        canvas.drawCentredString(0, 0, "SMITANGSHU BIO LAB")
        canvas.restoreState()

    def generate_pdf():
        doc = SimpleDocTemplate("report.pdf", pagesize=A4)
        styles = getSampleStyleSheet()

        elements = []

        elements.append(Paragraph("DNA Analysis Report", styles["Title"]))
        elements.append(Spacer(1, 10))

        table = Table([
            ["Length", len(dna)],
            ["GC %", round(gc, 2)],
            ["AT %", round(at, 2)]
        ])
        table.setStyle(TableStyle([
            ("GRID", (0,0), (-1,-1), 1, colors.black),
            ("BACKGROUND", (0,0), (-1,0), colors.lightgrey)
        ]))

        elements.append(table)
        elements.append(Spacer(1, 10))

        elements.append(Paragraph(f"<b>DNA:</b> {dna}", styles["Normal"]))
        elements.append(Paragraph(f"<b>Protein:</b> {protein}", styles["Normal"]))

        doc.build(elements,
                  onFirstPage=lambda c, d: [draw_border(c,d), draw_watermark(c,d)],
                  onLaterPages=lambda c, d: [draw_border(c,d), draw_watermark(c,d)])

        with open("report.pdf", "rb") as f:
            return f.read()

    pdf = generate_pdf()

    st.download_button("📄 Download Professional Report", pdf, "dna_report.pdf")


# ---------------- FOOTER ----------------
st.markdown("---")
st.markdown("💖 Made with 💕 by Smitangshu")