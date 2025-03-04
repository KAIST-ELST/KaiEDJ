#!/bin/bash
#$ -S /bin/bash
#$ -q intel2
#$ -pe intel2 64
#$ -l h=node11
#$ -v OMP_NUM_THREADS=64
#$ -cwd

echo "nodes list"
cat ${PE_HOSTFILE}
source /home/users1/hjun/.bashrc
source /home/users1/hjun/bashintel22_oneapi_py2.7_edmftf

pwd
echo "--------------------------------"
echo "jobid:    ${JOB_ID}"
echo "jname:    ${JOB_NAME}"
echo "nodes:    ${NNODES}"
echo "cores:    ${NSLOTS}"
echo "Nodefile: ${PE_HOSTFILE}"

echo "# time-start "`date`
echo "--------------------------------"

#node_list=./machines
#echo "machine-file : " $node_list
#echo "$ mpirun -np $NSLOTS hostname > $node_list"
#mpirun -np $NSLOTS hostname 
#cat $node_list
#julia --machine-file $node_list Jx_DMFT.jl -T fe_J_wannier.toml |& tee screen.out

julia Jx_DMFT.jl -T fe_J_wannier.toml 2>&1 | tee out

echo "------------------"
echo "# time-end " `date`
