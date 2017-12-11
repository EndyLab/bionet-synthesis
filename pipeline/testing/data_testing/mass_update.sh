
# $2 = field
# $3 = new data type

jq "$2 = \"$3\"" $1 | sponge $1
