sed -n '/# Gene/,/# Gene Description/ p' $1 | head -n -1 > tmp.md
md_to_json tmp.md > ID/${1::-2}json
rm tmp.md
sed -n '/# Gene Description/,$p' $1 > ID/${1::-3}-cut.md

