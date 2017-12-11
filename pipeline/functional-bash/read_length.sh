
# 1 file given 

parameter=$(./json_read.sh $1 $2)
if [ "$parameter" = "$3" ]; then
       fasta_unprocessed=${1::-5}
       fasta=$(echo "$fasta_unprocessed.fasta")
       ./sequence_length.sh $fasta write
fi


