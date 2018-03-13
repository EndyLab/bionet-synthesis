echo "Part type?"
read part_type
echo "Prefix?"
read prefix
echo "Suffix?"
read suffix

insert=$(printf "\t\"$part_type\": {\n\t\t\"prefix\": \"$prefix\",\n\t\t\"suffix\":\"$suffix\"\n\t},")

# Add to prefix: GAAGACATGAGGTCTCT
# Add to suffix: TGAGACCGCTATGTCTTC

echo "$insert" >> map.json


