#rm json_read.txt
touch json_read.txt
./find_recursive.sh $1 .json json_read.sh $2 write
cat json_read.txt | grep -o $3 | wc -l 
rm json_read.txt
