#!/bin/bash
domx=$1    
domy=$2    
domz=$3    
sizex=$4
sizey=$5
sizez=$6      
res=$7        
codeTimeLim=$8
wallTime=$9    

np=$((domx*domy*domz))
nnodes=$(((np-1) / 24 + 1))

# load any modules needed to compile or in any other way run this script
#module load intel impi cuda
#ml intel/16.0.3 mkl impi

# compile and configure the code
cd ..
./configure --with-problem=stratCoolingCosAD --enable-shearing-box --enable-fargo --with-eos=adiabatic --enable-mpi --enable-resistivity --enable-sts
make clean
./configure --with-problem=stratCoolingCosAD --enable-shearing-box --enable-fargo --with-eos=adiabatic --enable-mpi --enable-resistivity --enable-sts
make all 
cd bin

# run summitScript.sh, which makes job.sh
rm job.sh
chmod +x summitScript.sh
# args: 		    name of job   name of input file   nNodes   np  QOS     HH:MM:SS   
./summitScript.sh   athena_$np    athinput.txt         $nnodes  $np normal  $wallTime   

# make the input file
rm athinput.txt
./makeAthinput.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19}

# submit the job to SLURM
chmod +x job.sh
./job.sh
#sbatch -A ucb-testing job.sh

# take a look at the queue to make sure it went through
#squeue -u dagl4841
