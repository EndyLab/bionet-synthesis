jq "$2 = \"$3\"" $1 | sponge $1
echo "Updated file $1 section $2 with $3"
