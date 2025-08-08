# cd ./ic-bootstrapping/src_sage
# chmod +x run.sh
# ./run $filepath.sage$
#clear
result=$(sage "$1") 
echo "$result"
rm "$1.py"