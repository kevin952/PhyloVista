# phylo/views.py

from django.shortcuts import render
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Align import MultipleSeqAlignment
from .forms import PhyloForm
import ast
import string
from io import StringIO, BytesIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import parse
from Bio.Align import MultipleSeqAlignment
import matplotlib.pyplot as plt
import base64

def home(request):
    return render(request, 'phylo/base.html', {'form': PhyloForm()})


def parse_input(input_data, input_type):
    if input_type == 'distance_matrix':
        return parse_distance_matrix(input_data)
    elif input_type == 'fasta_sequence':
        return parse_fasta_sequence(input_data)
    else:
        return None
    

def parse_distance_matrix(matrix_str):
    # Convert the input string to a list of lists using ast.literal_eval
    data = ast.literal_eval(matrix_str)

    # Check if the data is a valid square matrix
    num_rows = len(data)

    # Generate labels as letters for each column
    labels = list(string.ascii_uppercase[:num_rows])

    return DistanceMatrix(labels, matrix=data)

def parse_fasta_sequence(fasta_str):
    # Create a file-like object from the string using StringIO
    fasta_file = StringIO(fasta_str)

    # Parse the FASTA sequences and create SeqRecord objects
    records = [SeqRecord(Seq(seq.seq), id=seq.id, description=seq.description) for seq in parse(fasta_file, format='fasta')]

    # Create a MultipleSeqAlignment object
    alignment = MultipleSeqAlignment(records)

    return alignment

def phylo_tree(request):
    print("Reached here!")
    if request.method == 'POST':
        form = PhyloForm(request.POST)

        if form.is_valid():
            input_type = form.cleaned_data['input_type']
            input_data = form.cleaned_data['input_data']
            algorithm = form.cleaned_data['algorithm']

            # Process the input data based on the selected input type
            data = parse_input(input_data, input_type)

            if input_type == 'distance_matrix':
                tree  = construct_phylo_tree_distance(data, algorithm)
            
            else:
                tree = construct_phylo_tree_fasta(data, algorithm)

            if tree:
                # Plot the phylogenetic tree
                plt.figure(figsize=(10, 8))
                Phylo.draw(tree, do_show=False)

                # Convert the plot to a base64-encoded image
                buffer = BytesIO()
                plt.savefig(buffer, format='png')
                plt.close()
                img_str = base64.b64encode(buffer.getvalue()).decode()

                # Render the result with the tree image
                return render(request, 'phylo/result.html', {'tree_image': img_str})

    else:
        form = PhyloForm()

    return render(request, 'phylo/base.html', {'form': form})



def construct_phylo_tree_fasta(data, algorithm):
    constructor = DistanceTreeConstructor()
    calculator = DistanceCalculator('identity')

    if algorithm == 'upgma':
        dm = calculator.get_distance(data)
        tree = constructor.upgma(dm)
    elif algorithm == 'neighbour_joining':
        dm = calculator.get_distance(data)
        tree = constructor.nj(dm)
    else:
        # Handle unsupported algorithm
        tree = None

    return tree

def construct_phylo_tree_distance(data, algorithm):
    # Convert the matrix to a list of lists
    matrix_list = [list(row) for row in data.matrix]

    # Create DistanceMatrix directly
    distance_matrix = DistanceMatrix(data.names, matrix=matrix_list)

    # Calculate the distance using the DistanceCalculator
    calculator = DistanceCalculator()
    dm = calculator.get_distance(distance_matrix)

    # Construct the phylogenetic tree
    constructor = DistanceTreeConstructor()
    if algorithm == 'upgma':
        tree = constructor.upgma(dm)
    elif algorithm == 'neighbour_joining':
        tree = constructor.nj(dm)
    else:
        # Handle unsupported algorithm
        tree = None

    return tree