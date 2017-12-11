# 1 directory
# 2 section to read
# 3 parameter you are looking for


echo "0" > sequence-total.txt

./find_recursive.sh $1 .json ./read_length.sh $2 $3

awk '{ sum += $1 } END { print sum }' sequence-total.txt
rm sequence-total.txt

