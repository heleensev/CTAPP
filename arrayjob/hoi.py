import timebuddy

from time import time
from time import sleep

def monitor_mongo(start, rt):
    st = time()

    buffer = 120

    # format runtime parameter to seconds
    rt = timebuddy.timeformat(rt)
    maxtime = rt.to_seconds() - buffer

    print('maxtime start {}'.format(str(maxtime)))
    if maxtime < buffer:
        maxtime = buffer
    # return to parent program if maxtime get exceeded
    timepast = int(time()-st) + int(start)
    print('timepast start: {}'.format(str(timepast)))
    maxtime = 20
    while timepast < maxtime:
        timepast = int(time()-st) + int(start)
        sleep(1)
    print('timepast: {}'.format(str(timepast)))
    print('maxtime {}'.format(str(maxtime)))
    print(time() - st)

monitor_mongo(2, '00:00:20')

def waiter():


    start = time()

    waittime = time()+ 15
    while time() < waittime:
        sleep(1)

    print(time()-start)
#
# waiter()