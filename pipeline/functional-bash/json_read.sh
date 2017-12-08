json_read=$( echo | jq -r $2 $1)
if [ "$3" = "write" ]; then
	echo "$json_read" >> json_read.txt
else
	echo "$json_read"
fi
