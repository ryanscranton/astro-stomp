import numpy
import subprocess
from optparse import OptionParser
import python_correlator

"""
DLS Wrapper code for running python_correlator. Loads catalogs, gets relavent 
columns and loops through all subfields in the survey.
"""

parser = OptionParser()
parser.add_option("-f","--foreground",dest="foreground",default="",
                  action="store",type="str",
                  help="Foreground file to load")
parser.add_option("-b","--background",dest="background",default="",
                  action="store",type="str",help="Background File to load")
parser.add_option("-o","--output_tag",dest="output_tag",
                  default="Wtheta_test",action="store",type="str",
                  help="Output tag for file names")
parser.add_option("-a","--auto_corr",dest="auto_corr",
                  default="false",action="store",type="str",
                  help="Run Auto-Correlation")
parser.add_option("-n","--n_randoms",dest="n_randoms",default=1,
                  action="store",type="int",help="# Randoms")
parser.add_option("--z_min",dest="z_min",default=0.5,
                  action="store",type="float",help="Minimum Redshift")
parser.add_option("--z_max",dest="z_max",default=1.2,
                  action="store",type="float",help="Maximum Redshift")
parser.add_option("--mag_min",dest="mag_min",default=23,
                  action="store",type="float",help="Min LBG Magnitude")
parser.add_option("--mag_max",dest="mag_max",default=27,
                  action="store",type="float",help="Max LBG Magnitude")
parser.add_option("-w","--weight",dest="weight",default="",
                  action="store",type="str",help="")
(options, args) = parser.parse_args()

def alphaM1Fun(mag):
	if options.weight == "Bouwens":
		return -0.27+1.21783*10**10*numpy.exp(
                    -0.921034*mag) #Bouwens 2007
	elif options.weight == "Steidel":
		return -0.4+9.94461*10**9*numpy.exp(
                    -0.921034*mag) #Steidel 1999
	elif options.weight == "Sawicki":
		return -0.74+1.1956*10**10*numpy.exp(
                    -0.921034*mag) #Sawicki 2006
	else:
            return numpy.ones(len(mag))

sub = 1
alpha = 2
delta = 3
R = 9
z_b = 13
fore_x = 29
fore_y = 30
back_x = 34
back_y = 35

class DummyJob(object):
    
    def __init__(self):
        self.returncode = -1

    def poll(self):
        return self.returncode

field_list = ["F1","F2","F3","F4","F5"]
subfield_list = ["p11","p12","p13",
                 "p21","p22","p23",
                 "p31","p32","p33"]

fore_file = open(options.foreground)
fore = []
fore_subfield = []
for line in fore_file:
	if line[0] == "#":
		continue
	tmp_list = line.split()
	fore.append([float(tmp_list[alpha]),float(tmp_list[delta]),
		     float(tmp_list[R]),float(tmp_list[z_b]),
		     float(tmp_list[fore_x]),float(tmp_list[fore_y])])
	fore_subfield.append(tmp_list[sub])
fore = numpy.array(fore)
fore_subfield = numpy.array(fore_subfield)
fore_mask = numpy.logical_and(
	numpy.logical_and(
		numpy.logical_and(fore[:,4]>1000,fore[:,4]<9000),
		numpy.logical_and(fore[:,5]>1000,fore[:,5]<9000)),
	numpy.logical_and(
		numpy.logical_and(fore[:,2]>20.0,fore[:,2]<24.0),
		numpy.logical_and(fore[:,3]>options.z_min,
				  fore[:,3]<options.z_max)))
fore = fore[fore_mask]
fore_subfield = fore_subfield[fore_mask]
fore = numpy.array([fore[:,0],fore[:,1],
		    numpy.ones(len(fore))]).transpose()

back_file = open(options.background)
back = []
back_subfield = []
for line in back_file:
	if line[0] == "#":
		continue
	tmp_list = line.split()
	back.append([float(tmp_list[alpha]),float(tmp_list[delta]),
		     float(tmp_list[R]),float(tmp_list[z_b]),
		     float(tmp_list[back_x]),float(tmp_list[back_y])])
	back_subfield.append(tmp_list[sub])
back = numpy.array(back)
back_subfield = numpy.array(back_subfield)
back_mask = numpy.logical_and(
	numpy.logical_and(
		numpy.logical_and(back[:,4]>1000,back[:,4]<9000),
		numpy.logical_and(back[:,5]>1000,back[:,5]<9000)),
	numpy.logical_and(back[:,2]>options.mag_min,
			  back[:,2]<options.mag_max))
back = back[back_mask]
back_subfield = back_subfield[back_mask]
back = numpy.array([back[:,0],back[:,1],
		    alphaM1Fun(back[:,2])]).transpose()

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
				print "Running Correlation on",field,subfield
				fore_name = ('/sandbox/morrison/LBG/'+field+'/'+
					     field+subfield+'_fore_tmp.ascii')
				back_name = ('/sandbox/morrison/LBG/'+field+'/'+
					     field+subfield+'_back_tmp.ascii')
				rand_name = ('/sandbox/morrison/LBG/'+field+'/'+
					     field+subfield+'_randoms.fiat')
				numpy.savetxt(fore_name,
					      fore[
						fore_subfield==field+subfield])
				numpy.savetxt(back_name,
					      back[
						back_subfield==field+subfield])
				job_list[idx] = subprocess.Popen(
					'python python_correlator.py'
					' -f'+fore_name+' -b'+back_name+
					' -r'+rand_name+
					' -oWtheta_'+field+subfield+'_'+
					options.output_tag+
					' -a'+options.auto_corr+
					' -n'+str(options.n_randoms),
					shell=True)
				break
while (job_list[0].poll()==None or
       job_list[1].poll()==None or
       job_list[2].poll()==None or
       job_list[3].poll()==None):
    continue
subprocess.Popen('rm /sandbox/morrison/LBG/F*/F*p*_fore_tmp.ascii')
subprocess.Popen('rm /sandbox/morrison/LBG/F*/F*p*_back_tmp.ascii')
