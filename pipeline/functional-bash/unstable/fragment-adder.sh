cat dictionary.csv | while read line
do
    ID=$(echo $line | cut -d',' -f1)    # get location
    gene_name=$(echo $line | cut -d',' -f2)    # get name
    fragtotal=$(cat Fragments_total.csv | grep -o $gene_name | wc -l)
    location=$(echo "../../../data/$ID/$ID.json")
    jq ".location.fragments = {\"${ID}_1\": \"\"}" $location | sponge $location
    for i in $(seq 2 $fragtotal);
        do
            NAME=$(echo ${ID}_${i})
            jq ".location.fragments |= . + {\"$NAME\": \"\" }" $location | sponge $location
            echo "$1"
        done

done
