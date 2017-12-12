
# Default data to find ../../data/

if [ "$1" = "data" ]; then
	data=$(echo "../../data/")
else
	data=$(echo ".")
fi

if [ "$2" = "edit" ]; then
	if [ "$3" = "replace" ]; then
		filetype=$(echo ".json")
		sponge=$(echo "NA")
		echo "Note: This section should largely only be for backfilling. Use with caution."
		echo "Please enter the section to edit."
		read edit_section
		echo "Please enter the final json edits."
		read json
		./find_recursive.sh $data $filetype json_write.sh $edit_section $json $sponge
		echo "Completed."
                echo "Please check output. If there are ANY errors, STOP NOW. If everything looks ok, enter 'sponge' to continue."
                read sponge
		./find_recursive.sh $data $filetype json_write.sh $edit_section $json $sponge
	fi
	if [ "$3" = "findif-replace" ]; then
		filetype=$(echo ".json")
		sponge=$(echo "NA")
		echo "Please enter the section to search for."
		read search
		echo "Please enter the term to search for."
		read term
		echo "Please enter section to edit."
		read edit_section
		echo "Please enter final json edits."
		read json 
		echo "$json" > final.txt
		./find_recursive.sh $data $filetype if_write.sh $search $term $edit_section $sponge
		echo "Completed."
		echo "Please check output. If there are ANY errors, STOP NOW. If everything looks ok, enter 'sponge' to continue."
		read sponge
		./find_recursive.sh $data $filetype if_write.sh $search $term $edit_section $sponge
		rm final.txt
	fi
fi


