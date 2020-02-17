# Compare_sequences_with_Dotplot
Compare sequences against other and represent the comparison as a dotplot

# Usage
blast_comparison_dotplot.sh -q query.fasta -s subject.fasta [-m 90] [-l 90] [-k] [-h]
[] means argument is optional

# Options
	-q	query file, fasta format
	-s	database file, fasta format
	-m	Maximum number of hits to consider per pair of query/subjects. Default 100
	-l	min_hit_len Default 200
	-k	activate log mode. Keep the two subfolders created for the analysis instead of removing them 
	-h	display this help

# Author
Damien Richard 2019

