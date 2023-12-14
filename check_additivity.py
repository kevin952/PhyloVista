import streamlit as st
import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
import numpy as np

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

def main():
    st.title("Additivity Checker App")

    # Input box for sequences or distance matrix
    input_type = st.selectbox("Select input type:", ["Sequences", "Distance Matrix"])

    if input_type == "Sequences":
        sequences_input = st.text_area("Enter sequences (one per line):", height=200)

        # Display current sequences if not cleared
        if st.button("Clear Sequences"):
            sequences_input = ""
        else:
            st.subheader("Current Sequences:")
            st.text(sequences_input)

        # Check additivity button
        if st.button("Check Additivity"):
            # Parse sequences from input
            sequences = [sequence.strip() for sequence in sequences_input.splitlines()]

            # Check if there are at least 2 sequences
            if len(sequences) < 2:
                st.error("You need at least 2 sequences to check additivity.")
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

        # Check additivity button for distance matrix
        if st.button("Check Additivity"):
            # Parse distance matrix from input
            distance_matrix = parse_distance_matrix(distance_matrix_input)

            if distance_matrix is None:
                st.error("Invalid input format for distance matrix. Please enter a valid NumPy matrix.")
                return

            # Display distance matrix
            st.subheader("Distance Matrix:")
            st.table(distance_matrix)

            # Check additivity and display result
            is_additive_result = is_additive(distance_matrix)
            st.success(f"Additivity Check Result: {is_additive_result}")

if __name__ == "__main__":
    main()
