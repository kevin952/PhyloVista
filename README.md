# PhyloVista: Real-time Phylogenetic Tree Builder

PhyloVista is an interactive Python application for real-time construction and visualization of phylogenetic trees based on biological sequences. The application supports inputting sequences either manually or from a FASTA file, offers two tree construction algorithms (Neighbor-Joining and UPGMA), and allows users to check the additivity of the input sequences.

## Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Running the Application](#running-the-application)
  - [Sequence Input](#sequence-input)
  - [Additivity Check](#additivity-check)
  - [Phylogenetic Tree Visualization](#phylogenetic-tree-visualization)
- [Examples](#examples)
- [Additional Information](#additional-information)
- [Contributing](#contributing)
- [License](#license)

## Requirements

To run PhyloVista, you need the following dependencies:

- Python 3.x
- Streamlit
- BioPython
- Matplotlib
- Pillow (PIL)
- NumPy
- Plotly

You can install these dependencies using the following command:

```
pip install -r requirements.txt
```

## Installation

- Clone the repository to your local machine:
``` git clone https://github.com/your-username/phylovista.git```

- Navigate to the project directory:
```cd phylovista```

- Install the required dependencies:
```pip install -r requirements.txt```

## Usage
- Running the Application
- Execute the following command to run the PhyloVista application:
```streamlit run phylovista_app.py```

- This will start the Streamlit development server, and the application will be accessible in your web browser a [http://localhost:8501/](http://localhost:8501/)

## Sequence Input
### Add Sequence-by-sequence:
- Enter sequences manually in the A, C, T, G format, one sequence per line.
- Click "Create Tree" to construct the phylogenetic tree.

### Add from FASTA:
- Upload a FASTA file containing biological sequences.
- The application constructs the phylogenetic tree based on the provided sequences.

### Additivity Check
- After adding sequences, you can check the additivity of the input:
- Click the "Check Additivity" button.
- The application displays the distance matrix and the result of the additivity check.

### Phylogenetic Tree Visualization
- The constructed phylogenetic tree is displayed in the main content area.
- Click on the image to enlarge it.
- Download the phylogenetic tree image using the provided link.

## Examples
### Example 1: Adding Sequences

```
ACGT
CGTA
TACG
```

### Example 2: Uploading FASTA
Upload a FASTA file with the following content:
```
>Sequence_1
ACGT
>Sequence_2
CGTA
>Sequence_3
TACG
```

## Additional Information
- The application provides information about the four-point condition in phylogenetic tree building.
- It explains the conditions and uses of the four-point condition for assessing tree topology and ensuring quality control.