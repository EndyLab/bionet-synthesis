ID=$( echo | jq -r '.id' $1)
name=$( echo | jq -r '.gene_name' $1)

echo "$ID"
echo "$ID,$name" >> dictionary.csv
