import streamlit as st
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
from PIL import Image
import random
import re

# Function to create a MultipleSeqAlignment object from a list of sequence strings
def create_alignment(sequences):
    seq_records = [SeqRecord(Seq(seq), id=f"Seq{i+1}") for i, seq in enumerate(sequences)]
    alignment = MultipleSeqAlignment(seq_records)
    return alignment

# Function to update the tree with a new sequence
def update_tree(sequence_data, algorithm):
    if len(sequence_data) == 0:
        return ["Please add a sequence to create a tree", None]

    # Create a MultipleSeqAlignment object
    alignment = create_alignment(sequence_data)

    # Create a DistanceCalculator object
    calculator = DistanceCalculator('identity')

    # # Choose the tree construction algorithm
    # if algorithm == 'nj':
    #     constructor = DistanceTreeConstructor(calculator)
    # elif algorithm == 'upgma':
    #     constructor = DistanceTreeConstructor(calculator, method='upgma')
    # else:
    #     return ["Invalid algorithm choice", None]
    constructor = DistanceTreeConstructor(calculator)
    updated_tree = constructor.build_tree(alignment)

    return ["Tree created", updated_tree]

# Streamlit App
st.title("Phylogenetic Tree Builder")

# # Algorithm choice dropdown
# algorithm = st.sidebar.selectbox("Choose Algorithm", ["Neighbor-Joining (nj)", "UPGMA"])

# Sidebar
st.sidebar.header("Choose Operation")
operation = st.sidebar.radio("", ["Add Sequence-by-sequence", "Add from FASTA"])

# Main Content
sequence_data = []
tree_nj = None

# Main Content
if operation == "Add Sequence-by-sequence":
    
    sequences_input = st.text_area("Enter sequences (one per line):", height=200)
        
    # Important note for sequences
    st.sidebar.markdown("**Important note:**")
    st.sidebar.markdown("Add sequences in A, C, T, G format, one sequence on a new line.")
    st.sidebar.markdown("**Example:**")
    st.sidebar.code("ACGT\nCGTA\nTACG")



    # Display current sequences if not cleared
    if st.button("Clear Sequences"):
        sequences_input = ""

    else:
        st.subheader("Current Sequences:")
        st.text(sequences_input)
    
    # Check additivity button
    if st.button("Create Tree"):
        
            # Check if sequences are in the valid format
        if not re.search(r'^[ACGTacgt\n]+$', sequences_input):
            st.error("Invalid sequence format. Sequences must consist of only 'A', 'C', 'G', or 'T' characters, one sequence per line.")
            
        else:  
            # Parse sequences from input
            sequences = [sequence.strip() for sequence in sequences_input.splitlines()]

            _, tree_nj = update_tree(sequences,"algorithm.lower()")

            # Display Tree
            st.subheader("Phylogenetic Tree")
            st.text("Click on the image to enlarge.")

            # Save the tree as an image
            image_path = "phylo_tree.png"
            plt.figure(figsize=(8, 8))
            Phylo.draw(tree_nj, axes=plt.gca())
            plt.axis('off')
            plt.savefig(image_path, bbox_inches='tight', pad_inches=0)
            plt.close()

            # Display the image
            image = Image.open(image_path)
            st.image(image, caption="Phylogenetic Tree")

        
else:
    sequence_data = []
    tree_nj = None
    st.subheader("Add from FASTA")
    fasta_file = st.file_uploader("Upload a FASTA File", type=["fa", "fasta"])
    
    if fasta_file is not None:
        new_sequences = []
        with open(fasta_file.name, 'r') as fasta_content:
            for line in fasta_content:
                if line.startswith(">"):
                    new_sequences.append("")
                else:
                    new_sequences[-1] += line.strip()

        # sequence_data.extend(new_sequences)
        # _, tree_nj = update_tree(sequence_data, "Seq_from_fasta")
        _, tree_nj = update_tree(new_sequences, "algorithm.lower()")

        # Display Tree
        if tree_nj is not None:
            st.subheader("Phylogenetic Tree")
            st.text("Click on the image to enlarge.")

            # Save the tree as an image
            image_path = "phylo_tree.png"
            plt.figure(figsize=(8, 8))
            Phylo.draw(tree_nj, axes=plt.gca())
            plt.axis('off')
            plt.savefig(image_path, bbox_inches='tight', pad_inches=0)
            plt.close()

            # Display the image
            image = Image.open(image_path)
            st.image(image, caption="Phylogenetic Tree")



