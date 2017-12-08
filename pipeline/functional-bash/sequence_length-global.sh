echo "0" > sequence-total.txt
./find_recursive.sh $1 .fasta sequence_length.sh write
awk '{ sum += $1 } END { print sum }' sequence-total.txt
rm sequence-total.txt
