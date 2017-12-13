
COUNTER=000000
mkdir nbdata
cat $1 | while read line
do
	let COUNTER=COUNTER+1
        CDS=$(echo $line | cut -d',' -f1)    # get CDS name
        sequence=$(echo $line | cut -d',' -f2)    # get sequence
        CDS=$(echo $CDS | sed -e 's/^"//' -e 's/"$//')
        sequence=$(echo $sequence | sed -e 's/^"//' -e 's/"$//')
        printf -v ID "%06d" "$COUNTER"
	full_id=$(echo "BBF10K-NOBUILD_$ID")
	echo "$full_id"
	mkdir nbdata/$full_id
        jq ".id = \"$full_id\"" template.json > nbdata/$full_id/$full_id.json
	jq ".gene_name = \"$CDS\"" nbdata/$full_id/$full_id.json | sponge nbdata/$full_id/$full_id.json
	echo ">$CDS" > nbdata/$full_id/$full_id.fasta
	echo "$sequence" >> nbdata/$full_id/$full_id.fasta
done
