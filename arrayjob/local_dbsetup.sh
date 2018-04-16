#!/bin/bash

START=$SECONDS
RT=$1

echo RUNTIME: $RT
echo START: $START
echo $SHELL

[ -f /hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/monrun.log ] && rm /hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/monrun.log

MONIP="$(ifconfig compute | grep "inet " | awk -F' ' '{print $2}')"
echo "MONIP=$MONIP" > /hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf
echo "MONHOST=$HOSTNAME" >> /hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf

module load sevpy

# mongod --config /hpc/local/CentOS7/dhl_ec/software/sevpy/1/etc/mongod.conf &

ping="$(/hpc/local/CentOS7/dhl_ec/software/ctapp/application/JobHandler/mongo_handler.py monitor_mongo ${START} ${RT})"
echo ${START}
echo PING $ping
# database should be terminated at this point, if monitor script finishes before runtime has elapsed,
# in this case we make the steps to shut down cleanly instead of getting a sigkill
if [ "${ping}" == "ok" ];
then
    mongod --shutdown --config /hpc/local/CentOS7/dhl_ec/software/sevpy/1/etc/mongod.conf

    rm /hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf

    echo "early shutdown at $(date)" > /hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/monrun.log
else
    echo "normal shutdown at $(date)" > /hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/monrun.log
    echo "no db running"
fi