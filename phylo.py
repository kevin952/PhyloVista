from Bio import Phylo
from Bio.SeqIO import parse
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

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
    with open(file_path, 'r') as fasta_file:
        records = list(parse(fasta_file, 'fasta'))
        for record in records:
            sequence_data.append(str(record.seq))
    return sequence_data


# Example sequences as strings
sequence_data = ["ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"]

# Create an initial tree
message, tree_nj = update_tree(sequence_data, "Seq1")

# Add a new sequence and update the tree
new_sequence = "ATGGCCATTGTAATGGTCCGCTGAAAGGGTGCCCGATAG"
sequence_data.append(new_sequence)
message, tree_nj = update_tree(sequence_data, "Seq2")

# Add a new sequence and update the tree
new_sequence = "ATGGCCATTGTAATGGTCCGCTGAAAGGGTGCCCGTTAG"
sequence_data.append(new_sequence)
message, tree_nj = update_tree(sequence_data, "Seq3")

# Display the updated tree
Phylo.draw_ascii(tree_nj)

# Add sequences from a FASTA file and update the tree
fasta_file_path = "alignment_file.fa"
new_sequences = add_sequences_from_fasta(fasta_file_path)
print(len(sequence_data))
sequence_data.extend(new_sequences)
print(len(sequence_data))
message, tree_nj = update_tree(sequence_data, "Seq4")

# Display the updated tree
Phylo.draw_ascii(tree_nj)







