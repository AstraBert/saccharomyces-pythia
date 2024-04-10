c=0
for i in test/*.fsa
do
    ((c++))
    python3 test/generate_test.py -i $i -o test/test${c}.csv
done 