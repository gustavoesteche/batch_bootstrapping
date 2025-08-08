# a scripting for running every test file on the repo
for var1 in $(find . -name 'test*'); do
    echo $var1
    ./run.sh "$var1"
    echo
done