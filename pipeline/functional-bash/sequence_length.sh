sequence=$(./sequence_retriever.sh $1)
sequencelength=${#sequence}


if [ "$2" = "write" ]; then
	echo "$sequencelength" >> sequence-total.txt
else
	echo "$sequencelength"
fi
