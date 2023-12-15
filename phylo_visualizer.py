import streamlit as st
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
from PIL import Image
import os
import re
import numpy as np
import base64


def calculate_distance_matrix(sequences):
    # Convert sequences to SeqRecord objects
    seq_records = [SeqRecord(Seq(seq), id=f"seq_{i}") for i, seq in enumerate(sequences)]

    # Convert SeqRecord objects to MultipleSeqAlignment
    alignment = MultipleSeqAlignment(seq_records)

    # Calculate Jukes-Cantor distance matrix
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Ensure the matrix is symmetric and fill symmetric values
    num_seqs = len(distance_matrix.names)
    distance_matrix_array = np.zeros((num_seqs, num_seqs), dtype=float)

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            value = distance_matrix[i, j]
            distance_matrix_array[i, j] = value
            distance_matrix_array[j, i] = value

    return distance_matrix_array

def is_additive(matrix):
    n = len(matrix)

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for l in range(k + 1, n):
                    left = matrix[i, j] + matrix[k, l]
                    right1 = matrix[i, k] + matrix[j, l]
                    right2 = matrix[i, l] + matrix[j, k]

                    if left >= right1 or left >= right2 or right1 != right2:
                        return False

    return True

def parse_distance_matrix(input_text):
    try:
        # Replace square brackets and split values
        values = input_text.replace("[", "").replace("]", "").split(',')

        # Convert to float
        distance_matrix = np.array(list(map(float, values)))

        # Determine the size of the square matrix
        size = int(len(distance_matrix) ** 0.5)

        # Reshape to a square matrix
        distance_matrix = distance_matrix.reshape((size, size))

        return distance_matrix
    except ValueError:
        return None
    

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

    # Choose the tree construction algorithm
    if algorithm == 'neighbor-joining (nj)':
        constructor = DistanceTreeConstructor(calculator, method='nj')
    elif algorithm == 'upgma':
        constructor = DistanceTreeConstructor(calculator, method='upgma')

    else:
        return ["Invalid algorithm choice", None]
    
    updated_tree = constructor.build_tree(alignment)

    return ["Tree created", updated_tree]

# Function to create a downloadable link for a file
def get_binary_file_downloader_html(file_path, download_link_text):
    with open(file_path, 'rb') as f:
        data = f.read()
    b64 = base64.b64encode(data).decode()
    return f'<a href="data:application/octet-stream;base64,{b64}" download="{file_path}">{download_link_text}</a>'


# Streamlit App
st.title("PhyloVista: Real-time Phylogenetic Tree builder")

# Algorithm choice dropdown
algorithm = st.sidebar.selectbox("Choose Algorithm", ["Neighbor-Joining (nj)", "UPGMA"])

# Sidebar
st.sidebar.header("Choose Operation")
operation = st.sidebar.radio("", ["Add Sequence-by-sequence", "Add from FASTA", "Check Additivity"])

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
        
    
        # Parse sequences from input
        sequences = [sequence.strip() for sequence in sequences_input.splitlines()]
        _, tree_nj = update_tree(sequences,algorithm.lower())

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
        st.image(image, caption=f"Phylogenetic Tree - {algorithm.lower()}")

        # Download button
        st.markdown(get_binary_file_downloader_html(image_path, 'Phylo_Tree_Image.png'), unsafe_allow_html=True)

        
elif operation == "Add from FASTA":
    sequence_data = []
    tree_nj = None
    st.subheader("Add from FASTA")
    fasta_file = st.file_uploader("Upload a FASTA File", type=["fa", "fasta"])

    # Sidebar note and example for FASTA sequences
    st.sidebar.markdown("**Note for FASTA sequences:**")
    st.sidebar.markdown("Upload a FASTA file containing biological sequences.")
    st.sidebar.markdown("**Example FASTA Format:**")
    st.sidebar.code(">Sequence_1\nACGT\n>Sequence_2\nCGTA\n>Sequence_3\nTACG")
    
    if fasta_file is not None:
        file_details = {"FileName":fasta_file.name,"FileType":fasta_file.type}
        st.write(file_details)

        with open(os.path.join("fasta_directory", fasta_file.name),"wb") as f: 
            f.write(fasta_file.getbuffer())         
            st.success("Saved File")
    
    if fasta_file is not None:
        new_sequences = []
        with open(os.path.join("fasta_directory", fasta_file.name), 'r') as fasta_content:
            for line in fasta_content:
                if line.startswith(">"):
                    new_sequences.append("")
                else:
                    new_sequences[-1] += line.strip()

        # sequence_data.extend(new_sequences)
        # _, tree_nj = update_tree(sequence_data, "Seq_from_fasta")
        _, tree_nj = update_tree(new_sequences, algorithm.lower())

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

            # Download button
            st.markdown(get_binary_file_downloader_html(image_path, 'Download Phylo Tree Image'), unsafe_allow_html=True)


else:

    # Input box for sequences or distance matrix
    input_type = st.selectbox("Select input type:", ["Sequences", "Distance Matrix"])

    if input_type == "Sequences":
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
        if st.button("Check Additivity"):
            # Check if sequences are in the valid format
            if not re.search(r'^[ACGTacgt\n]+$', sequences_input):
                st.error("Invalid sequence format. Sequences must consist of only 'A', 'C', 'G', or 'T' characters, one sequence per line.")
            else:
                # Parse sequences from input
                sequences = [sequence.strip() for sequence in sequences_input.splitlines()]

                # Check if there are at least 2 sequences
                if len(sequences) < 2:
                    st.error("You need at least 2 sequences to check additivity.")
                else:
                    # Check if sequences are of the same length
                    if len(set(len(seq) for seq in sequences)) > 1:
                        st.error("Sequences must be of the same length.")
                    else:
                        # Calculate distance matrix
                        distance_matrix = calculate_distance_matrix(sequences)

                        # Display distance matrix
                        st.subheader("Distance Matrix:")
                        st.table(distance_matrix)

                        # Check additivity and display result
                        is_additive_result = is_additive(distance_matrix)
                        st.success(f"Additivity Check Result: {is_additive_result}")

    elif input_type == "Distance Matrix":
        # Allow users to input a distance matrix
        st.subheader("Enter Distance Matrix:")
        distance_matrix_input = st.text_area("Enter the distance matrix (comma-separated values):", height=200)

        # Important note for distance matrix
        st.sidebar.markdown("**Important note:**")
        st.sidebar.markdown("Enter the distance matrix as a comma-separated values list. Rows and columns should be separated by new lines.")
        st.sidebar.markdown("**Example:**")
        st.sidebar.code("[0, 6, 13, 14],\n[6, 0, 15, 16],\n[13, 15, 0, 5],\n[14, 16, 5, 0]")

        # Check additivity button for distance matrix
        if st.button("Check Additivity"):
            # Parse distance matrix from input
            distance_matrix = parse_distance_matrix(distance_matrix_input)

            if distance_matrix is None:
                st.error("Invalid input format for distance matrix. Please enter a valid NumPy matrix.")
                
            else:
                # Display distance matrix
                st.subheader("Distance Matrix:")
                st.table(distance_matrix)

                # Check additivity and display result
                is_additive_result = is_additive(distance_matrix)
                st.success(f"Additivity Check Result: {is_additive_result}")

    
    # Static information about the four-point condition and its uses
    st.markdown("## Four-Point Condition in Phylogenetic Tree Building")
    st.markdown("The four-point condition is a criterion used in phylogenetic tree building to assess the reliability of the inferred relationships between taxa based on genetic distances.")
    st.markdown("### Conditions:")
    st.markdown("1. **Given three taxa A, B, and C:** If the distance from A to B is less than or equal to the sum of the distances from A to C and C to B, the four-point condition is satisfied.")
    st.markdown("2. **The condition must hold for all combinations of taxa:** This means that for every set of four taxa A, B, C, and D, the four-point condition should be satisfied.")
    st.markdown("### Uses:")
    st.markdown("- **Assessing Tree Topology:** Violation of the four-point condition may indicate inconsistencies in the inferred tree topology.")
    st.markdown("- **Quality Control:** Checking additivity is a step in ensuring the reliability of the phylogenetic tree.")
