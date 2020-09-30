for n in 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200
do
    ../project2 -Ijt -e 1.E-3 -n $n >>rotations_for_varying_n.csv
done
