import streamlit as st
from Bio import Phylo
from Bio.SeqIO import parse
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
from PIL import Image
import io

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

# Function to add sequences from a FASTA file to the existing sequence_data
def add_sequences_from_fasta(file_path):
    new_sequences = []
    with open(file_path, 'r') as fasta_file:
        records = list(parse(fasta_file, 'fasta'))
        for record in records:
            new_sequences.append(str(record.seq))
    return new_sequences

# Function to build a tree directly from a FASTA file
def build_tree_from_fasta(file_path):
    sequence_data = add_sequences_from_fasta(file_path)
    message, tree_nj = update_tree(sequence_data, "Seq_from_fasta")
    return message, tree_nj

# Streamlit App
st.title("Phylogenetic Tree Builder")

# Sidebar
st.sidebar.header("Choose Operation")
operation = st.sidebar.radio("", ["Add Sequence", "Add from FASTA", "Add FASTA to String Data"])

# Main Content
sequence_data = []
tree_nj = None

if operation == "Add Sequence":
    new_sequence = st.text_area("Enter New Sequence", "")
    if st.button("Add Sequence"):
        sequence_data.append(new_sequence)
        _, tree_nj = update_tree(sequence_data, f"Seq{len(sequence_data)}")

elif operation == "Add from FASTA":
    fasta_file = st.file_uploader("Upload a FASTA File", type=["fa", "fasta"])
    if fasta_file is not None:
        new_sequences = add_sequences_from_fasta(fasta_file)
        sequence_data.extend(new_sequences)
        _, tree_nj = update_tree(sequence_data, f"Seq{len(sequence_data)}")

elif operation == "Add FASTA to String Data":
    fasta_file = st.file_uploader("Upload a FASTA File", type=["fa", "fasta"])
    if fasta_file is not None:
        _, tree_nj = build_tree_from_fasta(fasta_file)

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


