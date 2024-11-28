import streamlit as st
from Bio import SeqIO
from io import StringIO

def read_fasta(file_content):
    """Reads a FASTA file content and returns a dictionary with sequence IDs as keys and sequences as values."""
    sequences = {}
    fasta_io = StringIO(file_content)
    for record in SeqIO.parse(fasta_io, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def calculate_gc_content(sequence):
    """Calculates the GC content of a given DNA sequence."""
    g_count = sequence.upper().count('G')
    c_count = sequence.upper().count('C')
    gc_content = (g_count + c_count) / len(sequence) * 100
    return gc_content

def find_snvs(sequence1, sequence2):
    """Finds single nucleotide variants (SNVs) between two DNA sequences."""
    snvs = []
    for i, (base1, base2) in enumerate(zip(sequence1, sequence2)):
        if base1 != base2:
            snvs.append((i + 1, base1, base2))  # Position is 1-indexed
    return snvs

def analyze_genome(sequences):
    """Analyzes a dictionary of sequences and returns a summary of the analysis."""
    analysis_results = {}
    for seq_id, sequence in sequences.items():
        sequence_length = len(sequence)
        gc_content = calculate_gc_content(sequence)

        analysis_results[seq_id] = {
            'length': sequence_length,
            'gc_content': gc_content
        }

    return analysis_results

def main():
    st.title("Genomic Data Analysis Tool")

    st.header("Upload FASTA File")
    fasta_file = st.file_uploader("Choose a FASTA file", type=["fna", "fasta"])
    if fasta_file:
        file_content = fasta_file.read().decode("utf-8")
        sequences = read_fasta(file_content)
        if sequences:
            st.subheader("Genome Analysis Results")
            analysis_results = analyze_genome(sequences)
            for seq_id, results in analysis_results.items():
                st.markdown(f"*Sequence ID*: {seq_id}")
                st.markdown(f"*Length*: {results['length']} bases")
                st.markdown(f"*GC Content*: {results['gc_content']:.2f}%")
                st.markdown("---")
        else:
            st.error("No sequences found in the file.")

    st.header("Single Nucleotide Variants (SNVs) Finder")
    seq1 = st.text_area("Enter Sequence 1")
    seq2 = st.text_area("Enter Sequence 2")

    if st.button("Find SNVs"):
        if seq1 and seq2:
            snvs = find_snvs(seq1, seq2)
            if snvs:
                st.subheader("Single Nucleotide Variants (SNVs):")
                for snv in snvs:
                    st.markdown(f"*Position {snv[0]}*: {snv[1]} -> {snv[2]}")
            else:
                st.info("No variants found.")
        else:
            st.error("Please enter sequences in both fields.")

    st.header("Analyze Evolutionary History")
    evo_file1 = st.file_uploader("Choose the first FASTA file", type=["fna", "fasta"], key="evo_file1")
    evo_file2 = st.file_uploader("Choose the second FASTA file", type=["fna", "fasta"], key="evo_file2")

    if st.button("Analyze Evolutionary History"):
        if evo_file1 and evo_file2:
            file_content1 = evo_file1.read().decode("utf-8")
            file_content2 = evo_file2.read().decode("utf-8")
            sequences1 = read_fasta(file_content1)
            sequences2 = read_fasta(file_content2)
            if sequences1 and sequences2:
                sequence1 = next(iter(sequences1.values()))
                sequence2 = next(iter(sequences2.values()))

                insertions = deletions = substitutions = 0
                length_diff = abs(len(sequence1) - len(sequence2))

                # Normalize the length of sequences
                if len(sequence1) > len(sequence2):
                    sequence2 += '-' * length_diff
                    insertions = length_diff
                elif len(sequence2) > len(sequence1):
                    sequence1 += '-' * length_diff
                    deletions = length_diff

                for base1, base2 in zip(sequence1, sequence2):
                    if base1 != base2:
                        if base1 == '-':
                            insertions += 1
                        elif base2 == '-':
                            deletions += 1
                        else:
                            substitutions += 1

                st.subheader("Evolutionary History Analysis")
                st.markdown(f"*Insertions*: {insertions}")
                st.markdown(f"*Deletions*: {deletions}")
                st.markdown(f"*Substitutions*: {substitutions}")
            else:
                st.error("No sequences found in one or both files.")
        else:
            st.error("Please upload both FASTA files.")

if __name__ == "__main__":
    main()
