import numpy
import subprocess
import time

class DummyJob(object):
    
    def __init__(self):
        self.returncode = -1

    def poll(self):
        return self.returncode

field_list = ["F1","F2","F3","F4","F5"]
subfield_list = ["p11","p12","p13",
                 "p21","p22","p23",
                 "p31","p32","p33"]

n = 0
clock = -1
hold_clock = -1
job_list = [DummyJob(),DummyJob(),DummyJob(),DummyJob()]
for field in field_list:
    for subfield in subfield_list:
        while (job_list[0].poll()==None and
               job_list[1].poll()==None and
               job_list[2].poll()==None and
               job_list[3].poll()==None):
            continue
        for idx,job in enumerate(job_list):
            if job.returncode != None:
                job_list[idx] = subprocess.Popen(
                    'python generate_randoms_wraper.py '+
                    ' -s'+field+subfield+
                    ' -o/sandbox/morrison/LBG/'+
                    field+'/'+field+subfield+
                    '_randoms.fiat -n 560000',
                    shell=True)
                print "Running",field,subfield
                break
while (job_list[0].poll()==None or
       job_list[1].poll()==None or
       job_list[2].poll()==None or
       job_list[3].poll()==None):
    continue
