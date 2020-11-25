export OMP_NUM_THREADS=1
echo 1 thread:
time ../../source/project4 -R 2 .01 2.4 -l 10 -c 1E5 >/dev/null
echo

OMP_NUM_THREADS=2
echo 2 threads:
time ../../source/project4 -R 2 .01 2.4 -l 10 -c 1E5 >/dev/null
echo

OMP_NUM_THREADS=4
echo 4 threads:
time ../../source/project4 -R 2 .01 2.4 -l 10 -c 1E5 >/dev/null
echo

OMP_NUM_THREADS=8
echo 8 threads:
time ../../source/project4 -R 2 .01 2.4 -l 10 -c 1E5 >/dev/null
