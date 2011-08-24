import numpy
import subprocess
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s","--subfield",dest="subfield",default="",
                  action="store",type="str",
                  help="Field To Create Randoms on")
parser.add_option("-o","--output_file",dest="output_file",default="test",
                  action="store",type="str",help="name of output catalog")
parser.add_option("-n","--n_randoms",dest="n_randoms", default=10000,
                  action="store",type="int",help="# of Randoms to Create")
(options, args) = parser.parse_args()

steps = numpy.ceil(options.n_randoms/(2000.0))

for i in xrange(int(steps)):
    job = subprocess.Popen('python generate_randoms.py -s '+options.subfield+
                           ' -o /sandbox/morrison/LBG/'+options.subfield[:2]+
                           '/'+options.subfield+'_randoms-'+str(i)+'.fiat'+
                           ' -n 2000',
                           shell=True)
    job.communicate()

cat_string = "cat "
for i in xrange(int(steps)):
    cat_string += ('/sandbox/morrison/LBG/'+options.subfield[:2]+
                   '/'+options.subfield+'_randoms-'+str(i)+'.fiat ')
cat_string += '> '+options.output_file
job = subprocess.Popen(cat_string,shell=True)
job.communicate()

for i in xrange(int(steps)):
    subprocess.Popen('rm /sandbox/morrison/LBG/'+options.subfield[:2]+
                     '/'+options.subfield+'_randoms-'+str(i)+'.fiat',shell=True)
