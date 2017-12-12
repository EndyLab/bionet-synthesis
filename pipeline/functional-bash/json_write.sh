if [ "$4" = "sponge" ]; then
	jq "$2 = $3" $1 | sponge $1
else
	jq "$2 = $3" $1
fi
echo "Updated file $1 section $2 with $3"
