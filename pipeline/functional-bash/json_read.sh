json_read=$( echo | jq -r $2 $1)
echo "$json_read"
