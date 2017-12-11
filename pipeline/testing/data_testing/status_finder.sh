status_of_json=$( echo | jq -r $2 $1)
echo "$status_of_json"
echo "$status_of_json" >> total_status.txt
