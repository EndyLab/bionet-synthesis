touch total_status.txt
find data/*/ -name "*.json" -exec ./status_finder.sh {} $1 \;

clear
cat total_status.txt | grep -o $2 | wc -l

rm total_status.txt
