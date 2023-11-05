from Bio import motifs
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from IPython.display import Image

# Read the sequences from the FASTA file
sequences = []
for record in SeqIO.parse("sequences.fasta", "fasta"):
    sequences.append(record.seq)

# Create a motif object with the sequences
m = motifs.create(sequences)

# Generate and display the WebLogo
m.weblogo("my_motif_logo.png")
Image(filename="my_motif_logo.png")
