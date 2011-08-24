import numpy
import subprocess
from optparse import OptionParser

"""
DLS specific random object generator on the sphere. Uses the DS9 regions file
that defines the exclusion regions to imprint the catalog gemoetry on the
random points.
"""

parser = OptionParser()
parser.add_option("-s","--subfield",dest="subfield",default="",
                  action="store",type="str",
                  help="Subfield To Create Randoms on")
parser.add_option("-o","--output_file",dest="output_file",default="test.fiat",
                  action="store",type="str",help="name of output catalog")
parser.add_option("-n","--n_randoms",dest="n_randoms", default=1000,
                  action="store",type="int",help="# of Randoms to Create")
(options, args) = parser.parse_args()

command = ('xy2sky -d /mnt/DLS/Stacks/2006/'+options.subfield[:2]+'/'+
           options.subfield[2:]+
           '/R/dlsmake.fits 0.5 0.5 0.5 10000.5 10000.5 10000.5 1000.5 0.5')
tmp_list = []
pos = subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
for k in pos.stdout:
    tmp_list.append([numpy.float(k.split()[0]),numpy.float(k.split()[1])])
tmp_array = numpy.array(tmp_list)
sign = 1.0
if (numpy.min(tmp_array[:,1]) < 0 and numpy.max(tmp_array[:,1]) < 0):
    sign = -1.0
z_min = numpy.min(numpy.cos(tmp_array[:,1]*numpy.pi/180.0))
z_max = numpy.max(numpy.cos(tmp_array[:,1]*numpy.pi/180.0))
ra_min = numpy.min(tmp_array[:,0])
ra_max = numpy.max(tmp_array[:,0])

dec_random = sign*numpy.arccos(numpy.random.uniform(
        z_min,z_max,size=options.n_randoms))*180/numpy.pi
ra_random = numpy.random.uniform(ra_min, ra_max, size=options.n_randoms)
randoms = numpy.array([ra_random, dec_random]).transpose()

radec_string = ""
for random in randoms:
    radec_string += str(random[0]) + ' ' + str(random[1]) + ' '
job = subprocess.Popen('sky2xy -j /mnt/DLS/Stacks/2006/'+
                       options.subfield[:2]+'/'+options.subfield[2:]+'/R/'+
                       'dlsmake.fits '+radec_string+'> /sandbox/morrison/LBG/'+
                       options.subfield[:2]+'/'+options.subfield+
                       '_tmp_file.fiat',
                       shell=True)
job.communicate()

tmp_file2 = open('/sandbox/morrison/LBG/'+options.subfield[:2]+'/'+
                 options.subfield+'_tmp_file2.fiat','w')
tmp_file2.writelines('# fiat 1.0\n')
tmp_file2.writelines('# ttype1 = ra\n')
tmp_file2.writelines('# ttype2 = dec\n')
tmp_file2.writelines('# ttype3 = x\n')
tmp_file2.writelines('# ttype4 = y\n')
tmp_file = open('/sandbox/morrison/LBG/'+options.subfield[:2]+'/'+
                 options.subfield+'_tmp_file.fiat')
for line in tmp_file:
    tmp_list = line.split()
    if numpy.float(tmp_list[4]) < 0.5 or numpy.float(tmp_list[4]) > 10000.5:
        continue
    if numpy.float(tmp_list[5]) < 0.5 or numpy.float(tmp_list[5]) > 10000.5:
        continue
    tmp_file2.writelines(tmp_list[0]+' '+tmp_list[1]+' '+
                         tmp_list[4]+' '+tmp_list[5]+'\n')
tmp_file.close()
tmp_file2.close()

job = subprocess.Popen('regfilter '+'/sandbox/morrison/LBG/'+
                       options.subfield[:2]+'/'+options.subfield+
                       '_tmp_file2.fiat /mnt/DLS/Stacks/2006/Regions/'+
                       options.subfield[:2]+'/'+options.subfield[2:]+'/'+
                       options.subfield+'R.all.reg  > '+options.output_file,
                       shell=True)
job.communicate()
subprocess.Popen('rm /sandbox/morrison/LBG/'+options.subfield[:2]+'/'+
                 options.subfield+'_tmp_file.fiat',shell=True)
subprocess.Popen('rm /sandbox/morrison/LBG/'+options.subfield[:2]+'/'+
                 options.subfield+'_tmp_file2.fiat',shell=True)
