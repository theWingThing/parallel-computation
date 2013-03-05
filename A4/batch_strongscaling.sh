#!/bin/bash
# This script is intepreted by the Bourne Shell, sh
#
# Documentation for SGE is found in:
# http://docs.oracle.com/cd/E19279-01/820-3257-12/n1ge.html
#
# Tell SGE which shell to run the job script in rather than depending
# on SGE to try and figure it out.
#$ -S /bin/bash
#
# Export all my environment variables to the job
#$ -V
# Tun the job in the same directory from which you submitted it
#$ -cwd
#
# Give a name to the job
#$ -N apf-strong-scale
#
#$ -pe orte 8
# Specify a time limit for the job
#$ -l h_rt=00:02:00
#
# Join stdout and stderr so they are reported in job output file
#$ -j y
#
# Run on the debug queue; only one node may be used
# To use more than one node, specify the "normal" queue
#$ -q normal.q
#
# Specifies the circumstances under which mail is to be sent to the job owner
# defined by -M option. For example, options "bea" cause mail to be sent at the 
# begining, end, and at abort time (if it happens) of the job.
# Option "n" means no mail will be sent.
#$ -m aeb
#
# *** Change to the address you want the notification sent to, and
# *** REMOVE the blank between the # and the $
# $ -M kmw015@eng.ucsd.edu
#

# Change to the directory where the job was submitted from
cd $SGE_O_WORKDIR

echo
echo " *** Current working directory"
pwd
echo
echo " *** Compiler"
# Output which  compiler are we using and the environment
gcc -v
echo
echo " *** Environment"
printenv

echo

echo ">>> Job Starts"
date

# Commands go here

mpirun -np 1 ./apf -n 400 -i 1000 -x 1 -y 1
mpirun -np 1 ./apf -n 400 -i 1000 -x 1 -y 1 -k
mpirun -np 2 ./apf -n 400 -i 1000 -x 1 -y 2
mpirun -np 2 ./apf -n 400 -i 1000 -x 1 -y 2 -k
mpirun -np 2 ./apf -n 400 -i 1000 -x 2 -y 1
mpirun -np 2 ./apf -n 400 -i 1000 -x 2 -y 1 -k
mpirun -np 4 ./apf -n 400 -i 1000 -x 1 -y 4
mpirun -np 4 ./apf -n 400 -i 1000 -x 1 -y 4 -k
mpirun -np 4 ./apf -n 400 -i 1000 -x 2 -y 2 
mpirun -np 4 ./apf -n 400 -i 1000 -x 2 -y 2 -k
mpirun -np 4 ./apf -n 400 -i 1000 -x 4 -y 1
mpirun -np 4 ./apf -n 400 -i 1000 -x 4 -y 1 -k
mpirun -np 8 ./apf -n 400 -i 1000 -x 1 -y 8
mpirun -np 8 ./apf -n 400 -i 1000 -x 1 -y 8 -k
mpirun -np 8 ./apf -n 400 -i 1000 -x 2 -y 4
mpirun -np 8 ./apf -n 400 -i 1000 -x 2 -y 4 -k
mpirun -np 8 ./apf -n 400 -i 1000 -x 4 -y 2
mpirun -np 8 ./apf -n 400 -i 1000 -x 4 -y 2 -k
mpirun -np 8 ./apf -n 400 -i 1000 -x 8 -y 1
mpirun -np 8 ./apf -n 400 -i 1000 -x 8 -y 1 -k


date
echo ">>> Job Ends"
