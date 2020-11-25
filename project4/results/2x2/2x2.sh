../../source/project4 -l 2 -E 1 -c 1e7 -e 10 >2x2.csv &

for c in 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9
do
    ../../source/project4 -l 2 -R 1 1 1 -c $c -e 10
done

