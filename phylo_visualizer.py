from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
import matplotlib.pyplot as plt
from io import BytesIO
import base64
import streamlit as st

# Function to calculate and visualize the phylogenetic tree
def visualize_phylogenetic_tree(sequences):
    # Find the length of the longest sequence
    max_length = max(len(seq) for seq in sequences)

    # Pad shorter sequences with gaps and trim longer sequences
    seq_records = [
        SeqRecord(Seq(seq.ljust(max_length, ' ')), id=f"Sequence_{i + 1}")
        for i, seq in enumerate(sequences)
    ]

    # Create a temporary alignment with the padded sequences
    temp_alignment = AlignIO.MultipleSeqAlignment(seq_records)

    # Read the existing alignment from a FASTA file
    alignment = AlignIO.read('/Users/kevindsouza/Documents/test-proj/alignment_file.fa', 'fasta')

    # Extend the alignment with the padded sequences
    alignment.extend(temp_alignment)

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)

    # Construct the phylogenetic tree using UPGMA algorithm
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    # Draw the phylogenetic tree
    fig = plt.figure()
    Phylo.draw(tree, do_show=False)

    # Save the figure to a BytesIO object
    buf = BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)

    # Convert the BytesIO object to a base64-encoded string
    image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')

    # Display the image in Streamlit
    st.image(f"data:image/png;base64,{image_base64}")

# Streamlit application
def main():
    st.title("Phylogenetic Tree Visualizer")

    # User input for sequences
    user_input = st.text_input("Enter a sequence:")
    sequences = [user_input]

    # Display the phylogenetic tree
    visualize_phylogenetic_tree(sequences)

    # Button to add more sequences
    if st.button("Add Sequence"):
        additional_sequence = st.text_input("Enter additional sequence:")
        sequences.append(additional_sequence)
        visualize_phylogenetic_tree(sequences)

if __name__ == "__main__":
    main()
