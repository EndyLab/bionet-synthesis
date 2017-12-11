echo "ID,Name" > dictionary.csv
find data/*/ -name "*.json" -exec ./dictionary.sh {} \;
