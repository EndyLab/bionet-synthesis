reply=$(cat final.txt)
echo "The reply is $reply"
json_read=$(./json_read.sh $1 $2)
if [ "$json_read" = "$3" ]; then
	echo "Match found for search terms in $1"
	if [ "$5" = "sponge" ] ; then
		jq "$4 = $reply" $1 | sponge $1
		echo "Updated file $1 section $2 with $3"
	else
		echo "Term to update current -" && jq $4 $1
		#jq "$4 = $reply" $1
	fi
else
	echo "No match found in $1"
fi
