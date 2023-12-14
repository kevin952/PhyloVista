import streamlit as st
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
from PIL import Image
import random

# Function to create a MultipleSeqAlignment object from a list of sequence strings
def create_alignment(sequences):
    seq_records = [SeqRecord(Seq(seq), id=f"Seq{i+1}") for i, seq in enumerate(sequences)]
    alignment = MultipleSeqAlignment(seq_records)
    return alignment

# Function to update the tree with a new sequence
def update_tree(sequence_data, identifier):
    if len(sequence_data) == 0:
        return ["Please add a sequence to create a tree", None]

    # Create a MultipleSeqAlignment object
    alignment = create_alignment(sequence_data)

    # Create a DistanceCalculator object
    calculator = DistanceCalculator('identity')

    # Neighbor-joining
    constructor_nj = DistanceTreeConstructor(calculator)
    updated_tree = constructor_nj.build_tree(alignment)

    return ["Tree created", updated_tree]

# Streamlit App
st.title("Phylogenetic Tree Builder")

# Sidebar
st.sidebar.header("Choose Operation")
operation = st.sidebar.radio("", ["Add Sequence-by-sequence", "Add from FASTA"])

# Main Content
sequence_data = []
tree_nj = None

# Main Content
if operation == "Add Sequence-by-sequence":
    sequence_data = []
    tree_nj = None
    st.subheader("Add Sequences One by One")
    st.text("Enter a new sequence and click 'Add Sequence' to update the tree.")
    key = random.randint(1, 1000)
    while st.button("Add Sequence", key):
        new_sequence = st.text_area("Enter New Sequence", "")
        key = random.randint(1, 1000)
        if new_sequence:
            sequence_data.append(new_sequence)
            _, tree_nj = update_tree(sequence_data, f"Seq{len(sequence_data)}")

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
        _, tree_nj = update_tree(new_sequences, "Seq_from_fasta")

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



