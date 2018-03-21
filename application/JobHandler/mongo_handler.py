#!/hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python

from pymongo import MongoClient, errors
from time import sleep, time
from os import path, remove
import timebuddy
import subprocess
import logging
import fire

logging.basicConfig(filename='mongohandler.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger()

config_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'
setupjob = '/hpc/local/CentOS7/dhl_ec/software/ctapp/arrayjob/dbsetup.sh'
shutjob = '/hpc/local/CentOS7/dhl_ec/software/ctapp/arrayjob/shutdb.sh'

class MongoHandlerError(Exception):
    # always exit program if error occurs
    def __init__(self, message):
        super(MongoHandlerError, self).__init__(message)


class MongoHandler:

    def __init__(self):
        self.runtime = str
        self.vmem = str
        self.mongo = bool
        self.mon_vmem = str

    def __get_host_info(self):
        with open(config_path, 'r') as mon_conf:
            lines = mon_conf.readlines()
        self.monhost = lines[0].split('=')[1].strip()
        self.monip = lines[1].split('=')[1].strip()
        logger.info(self.monip)
        logger.info(self.monhost)

    def __ping(self):
        try:
            client = MongoClient(self.monip)
            result = client.admin.command('ping')
            if result.get('ok'):
                return True
        except errors.ServerSelectionTimeoutError:
            "expected connection error"
            return False

    def launch_mongo(self, mongo, rt, vmem):
        self.runtime = rt
        self.vmem = vmem

        if not mongo:
            return

        print('launching mongo')

        # remove mongo config file from previous run
        def __monconf_check():
            try:
                remove(config_path)
            except OSError:
                pass

        def __qsub_job(rt, vmem):
            # add 20 min to runtime for database run
            monrt = timebuddy.timeformat(rt).add_min(60)
            logger.info('db runtime is: {}'.format(monrt))
            print('db runtime is: {}'.format(monrt))
            print('db vmem is: {}'.format(vmem))
            logger.info('db vmem is {}'.format(vmem))
            subprocess.call('qsub -l h_rt={0} -l h_vmem={1} {2} {0}'
                            .format('01:45:00', '40G', setupjob), shell=True)

        def __fetch_node_info():
            wait = 0
            while not path.exists(config_path):
                sleep(5)
                wait += 5
                if wait > 400:
                    raise MongoHandlerError('waiting too long for mongo config')
            self.__get_host_info()

        def __mongo_alive():
            wait = 0
            while not self.__ping():
                self.__ping()
                wait += 1
                if wait > 10:
                    logger.info('tried pinging 6 times')
                    raise MongoHandlerError('waiting too long for mongo alive')

        __monconf_check()
        __qsub_job(self.runtime, self.vmem)
        __fetch_node_info()
        __mongo_alive()

    def check_lock(self):
        self.__get_host_info()
        client = MongoClient(self.monip)

        def locks():
            cur_op = client.admin.command('currentOp')
            locks = cur_op.get('inprog')[0].get('locks')
            if locks:
                return True
        wait = 0
        while locks():
            try:
                sleep(1)
                wait += 1
                if wait > 500:
                    raise MongoHandlerError('waiting too long for lock')
            except errors.ConnectionFailure:
                logger.error('There is not server while checking lock running at {}'
                             .format(timebuddy.time()))

    def monitor_mongo(self, start, rt):
        buffer = 120
        # start timer of mongo script
        timer = time()
        # format runtime parameter to seconds
        rt = timebuddy.timeformat(rt)
        maxtime = rt.to_seconds() - buffer

        logger.info('maxtime start {}'.format(str(maxtime)))

        # return to parent program if maxtime get exceeded
        timepast = int(time()-timer) + int(start)
        logger.info('timepast start: {}'.format(str(timepast)))
        while timepast < maxtime:
            timepast = int(time()-timer) + int(start)
        logger.info('timepast: {}'.format(str(timepast)))
        logger.info('maxtime {}'.format(str(maxtime)))
        self.__get_host_info()
        if self.__ping():
            print('ok')
        else:
            print('nope')

    def stop_mongo(self, mongo, jids):
        if not mongo:
            return
        self.__get_host_info()
        monhost = self.monhost
        subprocess.call('qsub -hold_jid {} -l hostname={} ../arrayjob/shut_db.sh'
                        .format(','.join(jids), monhost), shell=True)


if __name__ == '__main__':
    fire.Fire(MongoHandler)
