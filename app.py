import streamlit as st
from Bio.Seq import Seq
from Bio import Entrez

# ------------------ CONFIG ------------------
st.set_page_config(page_title="DNA Analyzer", layout="wide")

Entrez.email = "smitangshudas23@gmail.com"  # IMPORTANT

# ------------------ TITLE ------------------
st.title("🧬 DNA Sequence Analyzer")
st.caption("Professional Bioinformatics Tool")

# ------------------ INPUT ------------------
dna = st.text_area("Enter DNA Sequence")

# ------------------ VALIDATION ------------------
def is_valid_dna(seq):
    return all(base in "ATCG" for base in seq.upper())

# ------------------ ANALYSIS ------------------
if st.button("Analyze"):

    if not dna:
        st.warning("Please enter a DNA sequence")
    elif not is_valid_dna(dna):
        st.error("Invalid DNA sequence (Only A, T, C, G allowed)")
    else:
        dna = dna.upper()
        seq = Seq(dna)

        reverse_comp = seq.reverse_complement()
        rna = seq.transcribe()
        protein = seq.translate()

        # ---------- SEQUENCES ----------
        st.subheader("🧪 Sequences")
        st.code(dna)
        st.code(str(reverse_comp))
        st.code(str(rna))
        st.code(str(protein))

        # ---------- MOTIF ----------
        st.subheader("🔍 Motif Finder")
        motifs_input = st.text_input("Enter motifs (comma separated)", "ATG")
        motifs = [m.strip().upper() for m in motifs_input.split(",") if m.strip()]

        for motif in motifs:
            positions = [i for i in range(len(dna)) if dna.startswith(motif, i)]
            st.success(f"{motif} → Positions: {positions}")

        # ---------- NCBI SEQUENCE SEARCH ----------
        st.subheader("🌐 NCBI Sequence Search")

        if st.button("Search NCBI Database"):
            try:
                with st.spinner("Searching NCBI..."):

                    handle = Entrez.esearch(
                        db="nucleotide",
                        term=dna[:50],   # use partial sequence
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

                        if results.strip():
                            st.subheader("📄 Retrieved sequences")

                            st.code(results[:1000])  # cleaner preview

                            st.download_button(
                                label="⬇️ Download full results",
                                data=results,
                                file_name="ncbi_results.fasta",
                                mime="text/plain"
                            )
                        else:
                            st.warning("Empty result returned")

            except Exception as e:
                st.error(f"Error: {e}")

        # ---------- ORF ----------
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

        for idx, (frame, start, end, seq) in enumerate(orfs):
            with st.expander(f"ORF {idx+1} (Length {len(seq)})"):
                st.code(seq)

# ------------------ NCBI GENE SEARCH ------------------
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
                    term=f"{gene_query}[Gene] AND Homo sapiens[Organism]",
                    retmax=3
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

                    if results.strip():
                        st.subheader("📄 Gene Results")

                        # show only headers nicely
                        lines = results.split("\n")
                        for line in lines:
                            if line.startswith(">"):
                                st.markdown(f"**{line[1:]}**")

                        st.code(results[:1000])

                        st.download_button(
                            label="⬇️ Download full gene data",
                            data=results,
                            file_name="gene_results.fasta",
                            mime="text/plain"
                        )
                    else:
                        st.warning("Empty result returned")

        except Exception as e:
            st.error(f"Error: {e}")

# ------------------ FOOTER ------------------
st.markdown("💖 Made with 💕 by Smitangshu")