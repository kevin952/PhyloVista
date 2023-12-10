# phylo/forms.py

from django import forms

class PhyloForm(forms.Form):
    input_type = forms.ChoiceField(
        choices=[('distance_matrix', 'Distance Matrix'), ('fasta_sequence', 'FASTA Sequence')],
        widget=forms.RadioSelect,
    )
    input_data = forms.CharField(widget=forms.Textarea)
    algorithm = forms.ChoiceField(
        choices=[('upgma', 'UPGMA'), ('neighbour_joining', 'Neighbour Joining')],
        widget=forms.RadioSelect,
    )
