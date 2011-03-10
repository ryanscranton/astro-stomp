#!/usr/bin/env python

"""
"""
from termcolor import colored
import stomp,math,pyfits,sys,os,time,socket,numpy
import load_sample,my_progress_bar
reload(load_sample)
import load_ascii_sample


############################################################
#### User defined parameters
############################################################

#### Parameters
r_p_min_kpc = 35.
r_p_max_kpc = 1000.


#### Path definition for CITA machines
#seed_path = '/'
hostname_string = socket.gethostname()
print '----',hostname_string
if hostname_string.find('trinity') > -1:
    seed_path = '/sandbox/sschmidt/STOMP/BRICE/'
else:
    seed_path = '/sandbox/sschmidt/STOMP/BRICE/'


## output filename
output_file_name = seed_path +'STOMP_OUTPUT/output_photoqso_low_35_1Mpc.fit'

## Map file (can be skipped)
#data_path = seed_path+'raid-cita/menard/DATA/SDSS/Galaxy/photo/DR7/photoz/eta.physics.ucdavis.edu/DR7/'
#map_file = data_path + 'stripe_photoz.hmap_basic'

## File of (photometric) sample to probe:
X_file = seed_path + 'DATA/photomqso_05_15_new.dat'

## File for spectroscopic sample:
target_file = seed_path+'DATA/dr7qso.fit'
target_z_flag = 'z'

############################################################
############################################################




print "Current output file: %s " %output_file_name
print "New output? [y/n]"
choice = sys.stdin.readline()
choice
if choice[0] == 'y':
    print "new output filename:"
    output_file_name = sys.stdin.readline().replace('\n','')
print colored("Output filename used: %s",'cyan') %output_file_name


## here I load the data
#    print 'Loading map: %s' %map_file
#    my_map = stomp.Map(map_file)
#    print colored("Map area: %3.2f deg^2" %my_map.Area() , 'cyan' )

print "Load X data? [1=Yes]"
choice = sys.stdin.readline()
if choice[0] == '1':
#CHANGED TO ASCII LOAD!!!!
    X_map = load_ascii_sample.load_ascii_sample(X_file,z_min=0.1)#,n_max=1000)#,z_max=1.5)

print "Load target data? [1=Yes]"
choice = sys.stdin.readline()
if choice[0] == '1':
    target_sample, target_map = load_sample.load_sample(target_file,z_min=0.1) #,n_max=100)#00)


#print "Median redshift:",numpy.median(X_sample.field('z'))

##############################
## We create the tree
X_tree_map = stomp.IndexedTreeMap()
for X in X_map:
    X_tree_map.AddPoint(X)


##############################
## We find pairs
indices = stomp.IndexVector()
neighbor_list = []
separation_list = []
physical_separation_list = []
z_X_list = []
z_target_list = []
master_index_list = []
master_list = []
n_neighbor_list = []

print 'Finding pairs and computing physical separation...'
print 'r_p min,max = ',r_p_min_kpc,r_p_max_kpc,' [kpc]'

i_target = -1L
n_target = target_map.size()
my_progress = my_progress_bar.progressBar(0,n_target, 50)
for target in target_map[:]:

    my_progress.updateAmount_and_write_if_needed(i_target)
    i_target += 1

    z_target = target_sample[target.Index()].field(target_z_flag)
    dist_to_target_Mpc = stomp.Cosmology_AngularDiameterDistance( numpy.double(z_target) )
    theta_min = r_p_min_kpc / (1e3 * dist_to_target_Mpc) * stomp.RadToDeg
    theta_max = r_p_max_kpc / (1e3 * dist_to_target_Mpc) * stomp.RadToDeg
    angular_bin_temp = stomp.AngularBin(theta_min,theta_max)

    X_tree_map.FindPairs(target, angular_bin_temp, indices)

    n_neighbor_list.append( indices.size() )
    z_target_list.append(z_target)


print "found %i targets" % len(z_target_list)
print "Average number of neighbors: %f" %numpy.average(n_neighbor_list)

n_neighbor_col = pyfits.Column(name='n_neighbor', format="J", array=n_neighbor_list)
z_target_col = pyfits.Column(name="z_target", format="E", array=z_target_list)
cols = pyfits.ColDefs([z_target_col, n_neighbor_col])
my_output = pyfits.new_table(cols)

print "Writing to %s..." % output_file_name
if os.path.exists(output_file_name):
    print 'Creating new file'
    os.remove(output_file_name)
my_output.writeto(output_file_name)


print "\a\a\a"

