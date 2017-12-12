# ./findprint_json.sh /data/ .parameter

find $1 -name "*$2" -exec ./$3 {} $4 $5 $6 $7 \;
