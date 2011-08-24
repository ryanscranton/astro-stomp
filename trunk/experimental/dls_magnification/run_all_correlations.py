import subprocess

print "Running z0.91.1 m2324"
job = subprocess.Popen('python run_correlations.py '+
                      '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
                      '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
                      '-o m2324_z0.91.1 -n50 --z_min=0.9 --z_max=1.1 '+
                      '--mag_min=23 --mag_max=24',shell=True)
job.communicate()
print "Running z0.91.1 m2425"
job = subprocess.Popen('python run_correlations.py '+
                      '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
                      '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
                      '-o m2425_z0.91.1 -n50 --z_min=0.9 --z_max=1.1 '+
                      '--mag_min=24 --mag_max=25',shell=True)
job.communicate()

print "Running z0.91.1 m2527"
job = subprocess.Popen('python run_correlations.py '+
                      '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
                      '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
                      '-o m2527_z0.91.1 -n50 --z_min=0.9 --z_max=1.1 '+
                      '--mag_min=25 --mag_max=27',shell=True)
job.communicate()
#
#print "Running z0.51.2 m2327 weighted"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o m2327_r50 -n50 --z_min=0.5 --z_max=1.2 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki',shell=True)
#job.communicate()

#print "Running z0.50.7 m2327 weighted"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o z0.50.7_r50 -n50 --z_min=0.5 --z_max=0.7 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki',shell=True)
#job.communicate()
#
#print "Running z0.70.9 m2327 weighted"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o z0.70.9_r50 -n50 --z_min=0.7 --z_max=0.9 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki',shell=True)
#job.communicate()
#
#print "Running z0.91.2 m2327 weighted"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o z0.91.2_r50 -n50 --z_min=0.9 --z_max=1.2 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki',shell=True)
#job.communicate()
#
#print "Running z0.91.1 m2327 weighted"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o z0.91.1_r50 -n50 --z_min=0.9 --z_max=1.1 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki',shell=True)
#job.communicate()
#
#print "Running z0.51.1 m2327 weighted"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o z0.51.1_r50 -n50 --z_min=0.5 --z_max=1.1 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki',shell=True)
#job.communicate()

# print "Running z0.40.5"
# job = subprocess.Popen('python run_correlations.py '+
#                        '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                        '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                        '-o z0.40.5_r50 -n50 --z_min=0.4 --z_max=0.5 '+
#                        '--mag_min=23 --mag_max=27 -wSawicki',shell=True)
# job.communicate()

# print "Running z0.40.5"
# job = subprocess.Popen('python run_correlations.py '+
#                        '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                        '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                        '-o auto_z0.40.5_r50 -n50 --z_min=0.4 --z_max=0.5 '+
#                        '--mag_min=23 --mag_max=27 -wSawicki -atrue',shell=True)
# job.communicate()

#print "Running z0.51.2 Auto Corr"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o auto_z0.51.2_r50 -n50 --z_min=0.5 --z_max=1.2 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki -atrue',shell=True)
#job.communicate()
#
#print "Running z0.50.7 Auto Corr"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o auto_z0.50.7_r50 -n50 --z_min=0.5 --z_max=0.7 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki -atrue',shell=True)
#job.communicate()
#
#print "Running z0.70.9 Auto Corr"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o auto_z0.70.9_r50 -n50 --z_min=0.7 --z_max=0.9 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki -atrue',shell=True)
#job.communicate()

#print "Running z0.91.2 Auto Corr"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o auto_z0.91.2_r50 -n50 --z_min=0.9 --z_max=1.2 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki -atrue',shell=True)
#job.communicate()
#
#print "Running z0.91.1 Auto Corr"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o auto_z0.91.1_r50 -n50 --z_min=0.9 --z_max=1.1 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki -atrue',shell=True)
#job.communicate()

#print "Running z0.51.1 Auto Corr"
#job = subprocess.Popen('python run_correlations.py '+
#                       '-f /sandbox/morrison/LBG/DLSforeground.fiat '+
#                       '-b /sandbox/morrison/LBG/DLSbDropouts.fiat '+
#                       '-o auto_z0.51.1_r50 -n50 --z_min=0.5 --z_max=1.1 '+
#                       '--mag_min=23 --mag_max=27 -wSawicki -atrue',shell=True)
#job.communicate()
