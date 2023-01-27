!#/bin/bash

# Primero descargo el genoma de Dianthus de la base de datos del NCBI Genome assembly (ASM2309106v1)
datasets download genome accession GCA_023091065.1 --include gff3,rna,cds,protein,genome,seq-report --filename GCA_023091065.1.zip

# A continuación se descomprime la carpeta obtenida: obtenemos una carpeta con el genoma en formato fasta. 
unzip GCA_023091065.1.zip

