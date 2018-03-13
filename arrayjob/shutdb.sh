#!/bin/bash
#$ -l h_rt=00:10:00
#$ -l h_vmem=5G
#$ -S /bin/bash
#$ -cwd
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N shutdb

module load sevpy

/home/dhl_ec/hseverin/deploy/application/JobHandler/mongo_handler.py check_lock

mongod --shutdown --config /hpc/local/CentOS7/dhl_ec/software/sevpy/1/etc/mongod.conf

rm /hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf
