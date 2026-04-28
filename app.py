import streamlit as st
from Bio.Seq import Seq
from Bio import Entrez
import re

# Required for NCBI
Entrez.email = "smitangshudas23@gmail.com"

# ---------------- PAGE ----------------
st.set_page_config(page_title="DNA Analyzer", layout="wide")

st.title("🧬 DNA Sequence Analyzer")
st.caption("Professional Bioinformatics Tool")

# ---------------- INPUT ----------------
st.subheader("📥 Enter DNA Sequence")

dna_input = st.text_area("Paste DNA sequence (A, T, G, C only)")

# Validate DNA
def is_valid_dna(seq):
    return re.fullmatch(r"[ATGCatgc]+", seq) is not None

dna = ""

if st.button("Analyze"):
    if not dna_input:
        st.warning("Please enter a DNA sequence")
    elif not is_valid_dna(dna_input):
        st.error("❌ Invalid DNA sequence! Use only A, T, G, C")
    else:
        dna = dna_input.upper()
        st.success("✅ Valid DNA sequence")

        # ---------------- BASIC ANALYSIS ----------------
        seq = Seq(dna)
        reverse_comp = str(seq.reverse_complement())
        rna = str(seq.transcribe())
        protein = str(seq.translate())

        st.subheader("🧬 Sequences")
        st.code(dna)
        st.code(reverse_comp)
        st.code(rna)
        st.code(protein)

# ---------------- MOTIF ----------------
if dna:
    st.subheader("🔍 Motif Finder")

    motifs_input = st.text_input("Enter motifs (comma separated)", "ATG")
    motifs = [m.strip().upper() for m in motifs_input.split(",") if m.strip()]

    for motif in motifs:
        positions = [i for i in range(len(dna)) if dna.startswith(motif, i)]
        st.success(f"{motif} → Positions: {positions}")

# ---------------- NCBI SEARCH (DNA) ----------------
st.subheader("🌐 NCBI Search (from DNA)")

if st.button("Search NCBI using DNA"):
    if not dna or len(dna) < 10:
        st.warning("Please analyze a valid DNA sequence first")
    else:
        try:
            with st.spinner("Searching NCBI..."):
                handle = Entrez.esearch(
                    db="nucleotide",
                    term=f"{dna[:50]}[All Fields]",
                    retmax=3
                )
                record = Entrez.read(handle)
                handle.close()

                ids = record["IdList"]

                if not ids:
                    st.warning("No matches found")
                else:
                    st.success(f"Found {len(ids)} matches")

                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        id=",".join(ids),
                        rettype="fasta",
                        retmode="text"
                    )
                    results = fetch_handle.read()
                    fetch_handle.close()

                    st.subheader("📄 Retrieved sequences")
                    st.code(results)

        except Exception as e:
            st.error(f"Error: {e}")

# ---------------- NEW: GENE SEARCH ----------------
st.subheader("🔎 NCBI Gene Search")

gene_query = st.text_input("Enter gene name (e.g. BRCA1)")

if st.button("Search Gene in NCBI"):
    if not gene_query:
        st.warning("Enter a gene name")
    else:
        try:
            with st.spinner("Searching gene in NCBI..."):
                handle = Entrez.esearch(
                    db="nucleotide",
                    term=gene_query,
                    retmax=2
                )
                record = Entrez.read(handle)
                handle.close()

                ids = record["IdList"]

                if not ids:
                    st.warning("No gene matches found")
                else:
                    st.success(f"Found {len(ids)} matches")

                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        id=",".join(ids),
                        rettype="fasta",
                        retmode="text"
                    )
                    results = fetch_handle.read()
                    fetch_handle.close()

                    st.subheader("📄 Retrived sequences")

                    st.text_area("Preview (first 2000 chars)", results[:2000], height=300)

                    st.download_button(
                        label="⬇️Download full results",
                        data=results,
                        file_name="ncbi_results.fasta",
                        mime="text/plain"
                     )

        except Exception as e:
            st.error(f"Error: {e}")

# ---------------- ORF ----------------
if dna:
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

    results = find_orfs(dna)

    for idx, (frame, start, end, seq) in enumerate(results):
        with st.expander(f"ORF {idx+1} (Length {len(seq)})"):
            st.code(seq)

# ---------------- FOOTER ----------------
st.markdown("💖 Made with 💕 by Smitangshu")