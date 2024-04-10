echo "TEST RESULTS" > test/test_results.stats

for i in {1..10}
do
    python3 scripts/predict_test.py -i test/test${i}.csv >> test/test_results.stats
done