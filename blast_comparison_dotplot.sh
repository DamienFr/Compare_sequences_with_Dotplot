#!/bin/bash


usage() { printf "Compare sequences against other and represent the comparison as a dotplot.\n\nUsage: $0 -q query.fasta -s subject.fasta [-m 90] [-l 90] [-k] [-h]
[] means argument is optional
\navailable arguments:
\n\t-q\tquery file, fasta format
\t-s\tdatabase file, fasta format
\t-m\tMaximum number of hits to consider per pair of query/subjects. Default 100
\t-l\tmin_hit_len Default 200
\t-k\tactivate log mode. Keep the two subfolders created for the analysis instead of removing them 
\t-d\tDebug mode For debuging purposes ONLY. Provide your own blast table [-d your_own_filtered_blast_table.tab]
\t-h\tdisplay this help\n" 1>&2; exit 1; }

log_mode=false

while getopts q:s:m:l:d:kh? opts; do
   case ${opts} in
      q) file=${OPTARG} ;;
      s) database=${OPTARG} ;;
      m) max_nb_hit=${OPTARG} ;;
      l) min_hit_len=${OPTARG} ;;
      d) debug_mode=${OPTARG} ;;
      k) log_mode=true ;;
      h|\?) usage;  exit 0 ;;
   esac
done

[ -z "$file" ]  &&  usage
[ -z "$max_nb_hit" ]  && max_nb_hit=100  && echo "-m (max hit nb) not specified, using default value 100" 
[ -z "$min_hit_len" ]  && min_hit_len=200  && echo "-l (min hit length) not specified, using default value 200" 

RED='\E[31;40m'
GREEN='\E[32;40m'

command -v makeblastdb >/dev/null 2>&1 || { printf "${RED}I require makeblastdb but it's not installed.\n to install it type into the terminal\n \"sudo apt-get install makeblastdb\" \n ${GREEN}" >&2; exit 1; }
command -v blastn >/dev/null 2>&1 || { printf "${RED}I require blastn but it's not installed.\n to install it type into the terminal\n \"sudo apt-get install blastn\" \n ${GREEN}" >&2; exit 1; }


####################### debut du nouveau script
file_name=$( echo $file | perl -ne '$_ = (split /\//,$_)[-1]; print')
database_name=$( echo $database | perl -ne '$_ = (split /\//,$_)[-1]; print')

if ! blastn -version | grep -q "2.9"; then echo "Code developped for Blast+ v.2.9.0, might or might not work with other blast+ versions"; fi

if ( ! grep -q "," $file ) && ( ! grep -q "," $database ) && ( ! grep -q '|' $file ) && ( ! grep -q '|' $database ); then
makeblastdb -in ${database} -dbtype nucl -parse_seqids -out db -title "DB_"

blastn  -num_alignments 1000000 -query ${file} -out resultat_blast.tmp -task blastn -db db -outfmt '10 std qlen slen gaps'
echo "qseqid,sseqid,pcent_length_QUERY,pcent_length_SUBJECT,pident Percentage of identical matches,Alignment length,Number of mismatch,Number of gap openings,Start of alignment in query,End of alignment in query,Start of alignment in subject(ref),End of alignment in subject(ref),evalue,bitscore,Query_seq_len,Sbjct_seq_len, gaps" > res_${file_name}_on_${database_name}
perl -F',' -ane '$qlength = 100*$F[3]/$F[12]; $slength = 100*$F[3]/$F[13]; $F[2] = $qlength . "," . $slength . "," . $F[2]; $toprint = join(",",@F);print $toprint' resultat_blast.tmp >> res_${file_name}_on_${database_name}

# on garde que les XX premieres lignes par read
#max_nb_hit=100
#min_hit_len=200
# env max_nb_hit=$max_nb_hit min_hit_len=$min_hit_len perl -F',' -ane 'if($. == 1){print}else{$h{$F[1] . $F[0]}++ ;if( $h{$F[1] . $F[0]} <= $ENV{max_nb_hit} && $F[5] > $ENV{min_hit_len}){print}}' /home/bacterio/Copy/01.BIOINFO/26.LM180_nanopore/03.check_pcur_on_chrom/res_pilon.fasta_249733nucl_splitted1_on_pilon.fasta_5166189nucl_splitted4

env max_nb_hit=$max_nb_hit min_hit_len=$min_hit_len perl -F',' -ane 'if($. == 1){print}else{$h{$F[1] . $F[0]}++ ;if( $h{$F[1] . $F[0]} <= $ENV{max_nb_hit} && $F[5] > $ENV{min_hit_len}){print}}'  res_${file_name}_on_${database_name} > res_${file_name}_on_${database_name}_filtered

if [ ${#debug_mode} -gt 0 ] ; then 
printf "\n\n\n######### ######### #########\nDebug mode activated, program behavior altered. Development purposes only\n######### ######### #########\n\n"
cp ${debug_mode} res_${file_name}_on_${database_name}_filtered
fi

nb_results=$(wc -l < res_${file_name}_on_${database_name}_filtered)

if [ $nb_results -gt 1 ] ;then

env file=$file_name database=$database_name perl -F',' -ane 'BEGIN{print "$ENV{database} $ENV{file}" . "\nNUCMER\n" };
if($.>1){ print ">" . $F[0] . " $F[1] $F[14] $F[15]\n$F[8] $F[9] $F[10] $F[11] 0 0 0\n0\n" }' res_${file_name}_on_${database_name}_filtered > ${file_name}_on_${database_name}.delta

mummerplot -size large -Q ${database} -R ${file} -t png -p ${file_name}_on_${database_name} ${file_name}_on_${database_name}.delta
# mummerplot -size large -Q ${database} -R ${file} -t postscript -p ${file_name}_on_${database_name} ${file_name}_on_${database_name}.delta
printf "\n#####################\nEnd of the pipeline\n#####################\n\nYour output files is\n${file_name}_on_${database_name}.png\n\n"

else 
printf "\n#####################\nNo lines in Blast output. Ending now ...\n#####################\n\n\n\n"
fi

else
printf "\nCommas or pipes in input file(s), cannot continue because scripts are parsing on commas and pipes makes mummerplot crash...\n"
fi

if ! $log_mode  ; then 
rm -f ${file_name}_on_${database_name}.gp
rm -f ${file_name}_on_${database_name}.delta
rm -f ${file_name}_on_${database_name}.fplot
rm -f ${file_name}_on_${database_name}.rplot
rm -f db*
rm -f resultat_blast.tmp
rm -f res_${file_name}_on_${database_name}
fi




