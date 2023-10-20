First, the reference genome of *Dianthus caryophyllus* was downloaded
using the reference ASM2309106v1.

``` bash
## First I download the Dianthus genome from the NCBI Genome assembly database (ASM2309106v1).
 
datasets download genome accession GCA_023091065.1 --include gff3,rna,cds,protein,genome,seq-report --filename GCA_023091065.1.$

# The folder obtained is then unzipped: we obtain a folder with the genome in fasta format.
unzip GCA_023091065.1.zip
```

## 1. Repeat annotation

### 1.1 **De novo** Repeat identification

We execute RepeatModeler.

``` bash
# Build Repeat Modeler database
ASM=/scratch/botany/jpicazo/annotation/repeatmodeler/dcar_asm.fasta
BuildDatabase -name dcaryophyllus $ASM

#Run RepeatModeler
RepeatModeler -database dcaryophyllus -pa 20 -LTRStruct >& run.out
```

### 1.2 Classification of unknown repetitive elements

``` bash
#  We add a prefix to each element found using seqkit:

cat dcaryophyllus-families.fa | seqkit fx2tab | awk '{ print "diaCarSQ01_"$0 }' | seqkit tab2fx > diaCarSQ01-families-prefix.fa
ls
grep -c ">" diaCarSQ01-families-prefix.fa
cat diaCarSQ01-families-prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > diaCarSQ01-families-prefix.fa.known
ls
cat diaCarSQ01-families-prefix.fa | seqkit fx2tab | grep  "Unknown" | seqkit tab2fx > diaCarSQ01-families-prefix.fa.unknown
less diaCarSQ01-families-prefix.fa
ls
grep -c ">" *.fa*
```

### 1.3 RepeatMasker

``` bash

# We make the directories were the outputs will be kept
mkdir -p logs 01_simple_out 02_eudicotyledon_out 03_known_out 04_unknown_out

# We Run repeatmasker using as input the reference genome. 
# -e ncbi: using NCBI BLAST as the search engine for identifying similar elements
# -a: The repeat alignments are writen as *.align* ouput file.
# -noint and -xsmall: for annotating simple repeats. 

RepeatMasker -pa 5 -a -e ncbi -dir 01_simple_out -noint -xsmall  diaCarSQ_edited.fa  2>&1 | tee logs/01_simplemask.log
```

Output:

1.  reference-genome.fasta.align: resulting repeated elements.

2.  reference-genome.fasta.cat.gz: complete RepeatMasker results.

3.  reference-genome.fasta.masked: the masked reference genome. I.e. a
    copy of the reference genome with repetitive elements for this
    round.

4.  reference-genome.fasta.out: tabbed list of annotated regions for
    this round of annotation/masking.

5.  reference-genome.fasta.tbl: summary of the overall genome
    composition based on repeat clusters.

The elements obtained are renamed:

``` bash
rename fasta simple_mask 01_simple_out/diaCarSQ*
rename .masked .masked.fasta 01_simple_out/diaCarSQ*
```

``` bash
#### Round 2: We annotate/mask eudicotyledons elements sourced from Repbase using the output from the 1st round of RepeatMasker
RepeatMasker -pa 16 -a -e ncbi -dir 02_eudicotyledons_out -nolow \
-species eudicotyledons 01_simple_out/diaCarSQ.edited.fa.masked.fasta 2>&1 | tee logs/02_eudicotyledonsmask.log

# Round 2: rename outputs
rename .fasta simple_mask 02_simple_out/diaCarSQ*
rename .masked .masked.fasta 02_simple_out/diaCarSQ*

#### round 3: annotate/mask known elements sourced from species-specific de novo repeat library using output froom 2nd round of RepeatMasker.
# For it -nolow \ -lib and the output of known families. 

RepeatMasker -pa 16 -a -e ncbi -dir 03_known_out -nolow \
-lib diaCarSQ_families.known \
02_eudicotyledons_out/diaCarSQ.eudicotyledons_mask.masked.fasta 2>&1 | tee logs/03_knownmask.log

# Round 3: rename outputs
rename eudicotyledons_mask.masked.fasta known_mask 03_known_out/diaCarSQ*
rename .masked .masked.fasta 03_known_out/diaCarSQ*


#### Round 4: We annotate/mask unknown elements sourced from species-specific de novo repeat library using output froom 3nd round of RepeatMasker
RepeatMasker -pa 16 -a -e ncbi -dir 04_unknown_out -nolow \
-lib diaCarSQ_families.unknown \
03_known_out/diaCarSQ.known_mask.masked.fasta 2>&1 | tee logs/04_unknownmask.log

# round 4: rename outputs
rename known_mask.masked.fasta unknown_mask 04_unknown_out/diaCarSQ*
rename .masked .masked.fasta 04_unknown_out/diaCarSQ*
```

It is necessary to make an out as a summary of everything obtained, a
complete summary of the repeated elements of the reference genome.

Round 1 of repeatMasker

``` bash
#!/bin/bash

# Variable genoma:
GENOMA=/mnt/nfs_storage9tb/export/data/rachel/dca_genome/GCA_023091065.1/diaCarSQ_edited.fa

# Repeatmasker:
RepeatMasker -pa 8 -a -e ncbi -dir 01_simple_out -noint -xsmall $GENOMA  2>&1 | tee logs/01_simplemask.log

# We rename the end of the fasta file by simple mask in the folder.
rename 's/fasta/simple_mask/' 01_simple_out/diaCarSQ_edited*
rename 's/\.masked/\.masked\.fasta/' 01_simple_out/diaCarSQ_edited*
```

Ronda 2 Repeat masker

``` bash
#!/bin/bash

########
# Variable genoma:
GENOMA=/mnt/nfs_storage9tb/export/data/rachel/dca_genome/GCA_023091065.1/diaCarSQ_edited.fa

# Variable of the simple repeat masker solution:
REPEAT01=/mnt/nfs_storage9tb/export/data/rachel/repeatmasker/01_simple_out/diaCarSQ_edited.fa.masked.fasta

# Repeatmasker:
RepeatMasker -pa 5 -a -e ncbi -dir 02_eudicotyledon_out -nolow -species eudicotyledons $REPEAT01
2>&1 | tee logs/02_eudycotiledonsmask.log

# We rename the end of the fasta file to eudcyt in the folder:
rename 's/fa/eudicot_mask/' 02_eudicotyledon_out/diaCarSQ_edited*
rename 's/\.masked/\.masked\.fasta/' 02_eudicotyledon_out/diaCarSQ_edited*
```

Round 3 Repeatmasker with known elements

``` bash
##### We check the masked.fasta output name
# Define the above output variable:
REPEAT02=/mnt/nfs_storage9tb/export/data/rachel/repeatmasker/02_eudicotyledon_out/diaCarSQ_edited.fa.masked.fasta.fasta.masked

# We run repeatMasker
RepeatMasker -pa 5 -a -e ncbi -dir 03_known_out -nolow -lib diaCarSQ_families.known $REPEAT02 2>&1 | tee logs/03_knownmask.log

# Rename the outputs:
rename 's/eudicots_mask/known_mask/' 03_known_out/diaCarSQ*
rename 's/\.masked/\.masked\.fasta/' 03_known_out/diaCarSQ*
```

Round 4 of Repeatmasker with the library of unknown elements:

``` bash
#!/bin/bash

# We select the masked genome from the previous round: 

REPEAT03=/mnt/nfs_storage9tb/export/data/rachel/repeatmasker/03_known_out/diaCarSQ_edited.known_mask.masked.fasta


# We run repeatmasker again.

RepeatMasker -pa 5 -a -e ncbi -dir 04_unknown_out -nolow -lib diaCarSQ_families.unknown $REPEAT03 2>&1 | tee logs/04_unknownmask.log

# We rename the output
rename 's/known_mask/unknown_mask.masked.fasta/' 04_unknown_out/diaCarSQ*
```

Round 5 to combine all results into a single file

``` bash
# Combine files with results - .cat.gz
cat 01_simple_out/diaCarSQ_edited.simple_mask.cat.gz \ 02_eudicotyledons_out/diaCarSQ_edited.eudicots_mask.cat.gz \
03_known_out/diaCarSQ_edited.known_mask.cat.gz \
04_unknown_out/diaCarSQ_edited.unknown_mask.cat.gz \
> 05_full_out/diaCarSQ_edited.full_mask.cat.gz

# Tabular files with all repetitions - .out
cat 01_simple_out/diaCarSQ_edited.simple_mask.out \
<(cat 02_eudicotyledons_out/diaCarSQ_edited.eudicots_mask.out | tail -n +4) \
<(cat 03_known_out/diaCarSQ_edited.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/diaCarSQ_edited.unknown_mask.out | tail -n +4) \
> 05_full_out/diaCarSQ_edited.full_mask.out

# We copy the tabular files for the single repeats - .out
cp 01_simple_out/diaCarSQ_edited.simple_mask.out 05_full_out/diaCarSQ_edited.simple_mask.out

# Combine tabular files for complex and interspaced repetitions - .out
cat 02_eudicotyledons_out/diaCarSQ_edited.eudicots_mask.out \
<(cat 03_known_out/diaCarSQ_edited.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/diaCarSQ_edited.unknown_mask.out | tail -n +4) \
> 05_full_out/diaCarSQ_edited.complex_mask.out

# Combine aligned repetitions for all repetitions - .align
cat 01_simple_out/diaCarSQ_edited.simple_mask.align \
02_eudicotyledons_out/diaCarSQ_edited.eudicots_mask.align \
03_known_out/diaCarSQ_edited.known_mask.align \
04_unknown_out/diaCarSQ_edited.unknown_mask.align \
> 05_full_out/diaCarSQ_edited.full_mask.align
```

Next, a summary table with all the results is generated:

``` bash
#!/bin/bash

# calculate the length of the genome sequence in the FASTA
allLen=`seqtk comp diaCarSQ_edited.fa | datamash sum 2`; #mi genoma de referencia

# calculate the length of the N sequence in the FASTA
nLen=`seqtk comp diaCarSQ_edited.fa | datamash sum 9`;

# Tabulate repeats per subfamily with total bp and proportion of genome masked
cat 05_full_out/diaCarSQ_edited.full_mask.out | tail -n +4 | awk -v OFS="\t" '{ print $6, $7, $11 }' |
awk -F '[\t/]' -v OFS="\t" '{ if (NF == 3) print $3, "NA", $2 - $1 +1; else print $3, $4, $2 - $1 +1 }' |
datamash -sg 1,2 sum 3 | grep -v "\?" |

# We rename diaCarSQ_edited.complex__mask.out to check if the programme went well.
awk -v OFS="\t" -v genomeLen="${allLen}" '{ print $0, $3 / genomeLen }' > 05_full_out/diaCarSQ_edited.complex__mask.out
```

### 1.4 MASKED genome FASTA file

First it is necessary to create a GFF3 file, the standard format for
genome annotation information. To do this, we will use the script
created by Daren Card [Daren’s
script](https://github.com/darencard/GenomeAnnotation/blob/master/rmOutToGFF3custom)

``` bash
#!/usr/bin/env bash
usage()
{
cat << EOF
rmOutToGFF3custom
Version 1.1 (2021-09-30)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or
desirable functioning.
This script converts the .out file from RepeatMasker to a GFF3 file. Note that the output
is probably not perfect GFF3, so beware with downstream applications. This script
emulates the rmOutToGFF3.pl script supplied with RepeatMasker but provides a fuller ID
("target=") for each element in column 9 of the GFF. This ID includes the matching element,
like rmOutToGFF3.pl, but also includes the repeat family: in the format <Family>/<Element>.
This change is because many matching elements produced from RepeatModeler have IDs that
provide no information about repeat family classification. Output is written to standard
output (SDOUT).
This script requires requires awk, which should be available on any standard Unix system.
rmOutToGFF3custom -o <RM.out> [-h] > <name.gff3>
OPTIONS:
        -h      usage information and help (this message)
    -o      RepeatMasker .out file
EOF
}

while getopts "ho:" OPTION
do
        case $OPTION in
                help)
                        usage
                        exit 1
                        ;;
        o)
            RMOUT=$OPTARG
            ;;
    esac
done

if [[ -z $RMOUT ]]
then
    usage
    exit 1
fi

cat <(echo "##gff-version 3") \
<(cat ${RMOUT} | tail -n +4 | \
awk -v OFS="\t" '{ if ($12 ~ /)/) print $5, "RepeatMasker", "dispersed_repeat", $6, $7, $1, $9, ".", "Target="$11"/"$10" "$14" "$13; \
else print $5, "RepeatMasker", "dispersed_repeat", $6, $7, $1, $9, ".", "Target="$11"/"$10" "$12" "$13 }' | \
awk -v OFS="\t" '{ if ($7 == "C") print $1, $2, $3, $4, $5, $6, "-", $8, $9; else print $0 }' | \
sort -k1,1 -k4,4n -k5,5n)
```

To execute the generated script we use the following orders:

``` bash

bash rmOutToGFF3custom.sh -o 05_full_out/diaCarSQ_edited.full_mask.out > 05_full_out/diaCarSQ.full_mask.gff3
bash rmOutToGFF3custom.sh -o 05_full_out/diaCarSQ_edited.simple_mask.out > 05_full_out/diaCarSQ.simple_mask.gff3
bash rmOutToGFF3custom.sh -o 05_full_out/diaCarSQ_edited.complex_mask.out > 05_full_out/diaCarSQ.complex_mask.gff3
```

Finally, the complete masked genome is generated with multiple simple
repeats masked as “soft” and complex interspaced repeats masked as
“hard”.

``` bash

#!/bin/bash
# We create a FASTA file of masks genome of soft single-repeat masks.
# -fi is the entry we want to give it. 
bedtools maskfasta -soft -fi diaCarSQ_edited.fa -bed 05_full_out/diaCarSQ.simple_mask.gff3 -fo 05_full_out/diaCarSQ.simple_mask.soft.fasta


# Hard-masked and complex-repeat genome
bedtools maskfasta -fi 05_full_out/diaCarSQ.simple_mask.soft.fasta -bed 05_full_out/diaCarSQ.complex_mask.gff3 -fo 05_full_out/diaCarSQ.simple_mask.soft.complex_mask.hard.fasta
# If they are capitalised, they are gene regions (a priori) or structural DNA. And in lower case there are single elements and N interspersed elements.  


### Masking the complexes

# A fully masked genome is now created with hard masking to check that the repeats have been masked.

bedtools maskfasta -fi diaCarSQ.fasta -bed 05_full_out/diaCarSQ.full_mask.gff3 \
-fo 05_full_out/diaCarSQ.full_mask.hard.fasta
```

# 2. LTR Harvest

Before running the LTR harvest it is necessary to generate an index with
the suffixerator program.

``` bash
# We first activate the conda ambiance  
conda activate ltr_harvester

# Caryophyllus
gt suffixerator -db diaCarSQ_edited.fa -indexname diaCarSQ_edited_index.fa -tis -suf -lcp -des -ssp -sds -dna

# Broteri
gt suffixerator -db diaBro01.fasta -indexname diaBro_index.fa -tis -suf -lcp -des -ssp -sds -dna

# In the folder the files with some of the terminations are not generated, but it is not bad, in the github it specifies that in the last versions it is not generated.
```

We then proceed to launch ltr harvest with the following parameters:

``` bash
# Caryophyllus
gt ltrharvest -index diaCarSQ_edited_index.fa -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 > genome.harvest.scn

# Broteri:
gt ltrharvest -index diaBro_index.fa -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 > genome.diaBro.harvest.scn
```

-similar: minimum similarity between two LTRs (up to 90, dafault 85).

-vic: number of nucleotides to be searched for TSD or other motif.

-seed: the minimum length for repeats. Those below the threshold are not
counted.

-   minlenltr: minimum length of the LTR.

-   maxlenltr: maximum length.

-   mintsd: a TSD (duplication target site) with this minimum length is
    searched for.

-   maxtsd: maximum TSD length.

-motif: specify two nucleotides from the start of the motif and two from
the end of each LTR.

# 3. LTR Retriever

With the result obtained in the ltr harvest we can proceed to run the
ltr retriever to obtain the intact elements.

``` bash
# Caryophyllus:
LTR_retriever -genome diaCarSQ_edited.fa  -inharvest genome.harvest.scn -threads 16 > retriever_SQ_01.txt
#  We maintain the default parameters using a mutation rate of 1.3 x 10 -8, that of Oryzia sativa.

# Broteri:
LTR_retriever -genome diaBro01.fasta -inharvest genome.diaBro.harvest.scn -threads 8 > retriever_bro_01.txt

```

# 4. TEsorter

``` bash
conda activate ltr_harvester


## Extract sequences from LTR_retriever in LTR_restreiever folder

## DIANTHUS BROTERI

LTR_retriever.py get_full_seqs diaBro01.fasta > ltr_intacDiaBro.fa



### TE clasiffiy

TEsorter ltr_intacDiaBro.fa -db rexdb-plant -p 22



## DIANTHUS CARYOPHYLLUS


LTR_retriever.py get_full_seqs diaCarSQ_edited.fa > ltr_intacDiaCarSQ.fa



### TE clasiffiy

TEsorter ltr_intacDiaCarSQ.fa -db rexdb-plant -p 22
```
