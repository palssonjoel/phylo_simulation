from Bio import SeqIO

# Convert PHYLIP alignment to FASTA
SeqIO.convert("Aswine_Arnsberg_6554_1979_merged.fasta.phy", "phylip", "alignment.fasta", "fasta")
