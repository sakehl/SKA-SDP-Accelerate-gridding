#!/bin/bash
#SBATCH -t 0:10:00
#SBATCH -N 1
#SBATCH -p gpu_short

NPROC=`nproc --all`
M=1
M2=1
date=$(date '+%Y_%m_%d')TEST

cd $HOME/lvandenhaak/SKA-SDP-Accelerate-gridding/
mkdir "$TMPDIR"/lvdh
cp -r data "$TMPDIR"/lvdh

# for N in 100 1000 5000
# do

#     for value in `seq 1 $M`
#     do
#         echo "Job number $value for CPU n=$N" |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
#         /usr/bin/time --verbose stack run -i "$TMPDIR"/lvdh/data -n $N |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
#     done

#     for value in `seq 1 $M`
#     do
#         echo "Job number $value for GPU n=$N" |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
#         /usr/bin/time --verbose stack run -g -i "$TMPDIR"/lvdh/data -n $N |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
#     done
# done


################ SEQUENCES
cd $HOME/lvandenhaak/SKA-SDP-Accelerate-griddingSEQ/
export LD_LIBRARY_PATH=/home/gkeller/llvm-4.0/lib/:$LD_LIBRARY_PATH
export PATH=/home/gkeller/llvm-4.0/bin/:$PATH

#GPU
for N in 10000
do
	for value in `seq 1 $M`
    do
        echo "Job number $value for GPU Sequences n=$N" |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
        /usr/bin/time --verbose stack run -- -chunks 30000 -i "$TMPDIR"/lvdh/data -n $N -g |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
    done

    for value in `seq 1 $M`
    do
        echo "Job number $value for GPU Foreign Sequences n=$N" |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
        /usr/bin/time --verbose stack run -- -chunks 30000 -i "$TMPDIR"/lvdh/data -n $N -g -for |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
    done
done

#CPU
for N in 10000
do
	for value in `seq 1 $M2`
    do
        echo "Job number $value for GPU Sequences n=$N" |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
        /usr/bin/time --verbose stack run -- -chunks 30000 -i "$TMPDIR"/lvdh/data -n $N |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
    done

    for value in `seq 1 $M2`
    do
        echo "Job number $value for GPU Foreign Sequences n=$N" |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
        /usr/bin/time --verbose stack run -- -chunks 30000 -i "$TMPDIR"/lvdh/data -n $N -for |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
    done
done

# #CPU foreign one extra
# for N in 10000
# do
#     for value in `seq 1 $M2`
#     do
#         echo "Job number $value for GPU Foreign Sequences n=$N" |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
#         /usr/bin/time --verbose stack run -- -chunks 30000 -i "$TMPDIR"/lvdh/data -n $N -for |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
#     done
# done

######################## PYTHON
# cd $HOME/lvandenhaak/crocodile
# . pipenv/bin/activate

# for N in 100 1000 5000 10000
# do
#     for value in `seq 1 $M`
#     do
#         echo "Job number $value for python n=$N parallel $NPROC" |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
#         /usr/bin/time --verbose python3 scripts/image_dataset.py --theta 0.008 --lambda 300000 --wkern "$TMPDIR"/lvdh/data/SKA1_Low_wkern2.h5 --akern "$TMPDIR"/lvdh/data/SKA1_Low_akern3.h5 --max $N --image "test" -N $NPROC "$TMPDIR"/lvdh/data/SKA1_Low_quick.h5 |& tee -a $HOME/lvandenhaak/$date\_results$N.txt
#     done
# done