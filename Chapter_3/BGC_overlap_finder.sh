#!/bin/bash

# This script takes in a list of query contigs (ensure '>' at start of line is not included) and returns contigs that overlap with these contigs
# The merging of all plate well assemblies should be done outside of this script since it is only required once. Code is below (commented out).
# Two parameters that need to be adjusted in this script is the path to the vector sequence and the path to the reference assembly. It is easier to hard-code these rather than specify them each time the script is run.
# The parameters are as follows:
# $1 is the list of query contigs
# $2 is the fasta file these contigs originate from
# $3 is the working directory the output is written to
# Run the script as below:
# BGC_overlap_finder.sh /path/to/query_contig_list.txt /path/to/assembly_with_query_contigs.fasta /path/to/working_directroy

# Create working directory
mkdir $3
WORKDIR=$3
cd $WORKDIR

# Extract contigs from specified file
seqkit grep -n -f $1 $2 > ./query_contigs.fasta

# Map vector sequence to query contigs. Make sure the correct digested vector sequence is specified
minimap2 -x map-pb -a /nfs/scratch/nowakvi/Metagenome_v_cosmid/pWEB_smaI.fasta ./query_contigs.fasta > ./vector_to_query.sam

# Trim vector sequences from query sequences
~/Get_soft_clip_alt_3_20220610.py -s $WORKDIR/vector_to_query.sam -o $WORKDIR/trimmed_query_contigs.fasta

# Extract first and last 1000 bases from all vector-trimmed contigs. Have to rename to get unique names for end-seqs (*_2 is tail end)
seqkit subseq -r 1:1000 ./trimmed_query_contigs.fasta > ./trimmed_query_contig_end_seqs.fasta
seqkit subseq -r -1000:-1 ./trimmed_query_contigs.fasta >> ./trimmed_query_contig_end_seqs.fasta
seqkit rename ./trimmed_query_contig_end_seqs.fasta > ./trimmed_query_contig_end_seqs_renamed.fasta

# Should merge all contig files and inlcude filename
# From https://www.unix.com/unix-for-dummies-questions-and-answers/242665-append-file-name-fasta-file-headers-linux.html
# awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1'
# The below command has to be run in the directory where all the fasta files are, otherwise file path is included
# for i in *.fasta; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' $i >> ./all_Sp_contigs_renamed.fasta; done

# Map end-seqs to merged plate well assemblies
minimap2 -x map-pb -t 8 -a /nfs/scratch/nowakvi/Metagenome_v_cosmid/Sp/all_scaffolds/all_Sp_contigs_renamed.fasta ./trimmed_query_contig_end_seqs_renamed.fasta > ./trimmed_query_contig_end_seqs_renamed_to_merged_plate_wells.sam

# Identify overlap contigs and extract from merged plate well assemblies
samtools view -b ./trimmed_query_contig_end_seqs_renamed_to_merged_plate_wells.sam > ./trimmed_query_contig_end_seqs_renamed_to_merged_plate_wells.bam
samtools sort ./trimmed_query_contig_end_seqs_renamed_to_merged_plate_wells.bam > ./trimmed_query_contig_end_seqs_renamed_to_merged_plate_wells_sorted.bam
samtools index ./trimmed_query_contig_end_seqs_renamed_to_merged_plate_wells_sorted.bam
samtools view -F 4 ./trimmed_query_contig_end_seqs_renamed_to_merged_plate_wells_sorted.bam | awk -F'\t' '{print $3}' > ./list_of_overlap_contigs.txt
seqkit grep -n -f ./list_of_overlap_contigs.txt /nfs/scratch/nowakvi/Metagenome_v_cosmid/Sp/all_scaffolds/all_Sp_contigs_renamed.fasta > ./overlap_contigs.fasta
