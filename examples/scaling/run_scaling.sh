#!/bin/bash
# MPI scaling benchmark driver.
#   strong scaling : fixed 6.4e7 photons, np = 1..64 (equal-share) / 2..64 (master-slave)
#   weak scaling   : 1e6 photons per rank
# Writes results.csv: test,mode,np,nphotons,wall_seconds
cd "$(dirname "$0")"
EXE=../../MoCafe.x
STRONG_N=6.4e7
WEAK_PER=1.0e6
OUT=results.csv
echo "test,mode,np,nphotons,wall_seconds" > $OUT

run_one () {  # $1=test $2=mode(ms|eq) $3=np $4=nphotons
    local ms=.false.; [ "$2" = ms ] && ms=.true.
    sed -e "s/par%no_photons      = .*/par%no_photons      = $4/" \
        -e "s/par%use_master_slave = .*/par%use_master_slave = $ms/" \
        scale_base.in > run.in
    local t0=$(date +%s.%N)
    mpirun -np $3 $EXE run.in > last.log 2>&1
    local t1=$(date +%s.%N)
    local dt=$(echo "$t1 - $t0" | bc)
    echo "$1,$2,$3,$4,$dt" >> $OUT
    echo "  $1 $2 np=$3 N=$4 : ${dt}s"
    rm -f scale*.h5
}

echo "=== strong scaling, equal-share ==="
for np in 1 2 4 8 16 32 64; do run_one strong eq $np $STRONG_N; done
echo "=== strong scaling, master-slave ==="
for np in 2 4 8 16 32 64; do run_one strong ms $np $STRONG_N; done
echo "=== weak scaling (1e6 photons per rank) ==="
for np in 1 2 4 8 16 32 64; do
    n=$(python3 -c "print(f'{$WEAK_PER*$np:.1e}')")
    run_one weak eq $np $n
done
for np in 2 4 8 16 32 64; do
    n=$(python3 -c "print(f'{$WEAK_PER*$np:.1e}')")
    run_one weak ms $np $n
done
rm -f run.in last.log
echo "done -> $OUT"
