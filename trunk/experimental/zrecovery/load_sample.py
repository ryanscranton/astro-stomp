def load_sample(sample_file,x='RA',y='DEC',n_max=10000000,z_min=0,z_max=99):
    
    """ load SDSS data from fits files and create a target_map and a galaxy_map
    """

    import stomp
    import math
    from termcolor import colored
    import pyfits
    import sys
    import os
    import time
    import numpy

    
    ## initialize the oututs:
    sample_map = 0
    print "Loading sample file: %s" %sample_file
    sample_in = pyfits.getdata(sample_file)
    n_object_found = sample_in.size
    print colored("Found: %i objects" % n_object_found , 'cyan')
    

    # #############################
    # load sample of a stomp map
    
    sample_map = stomp.IAngularVector() #sample_map = stomp.WAngularVector()
    n_sample_read = -1
    n_sample_loaded = 0
    t0 = time.time()
    for sample in sample_in[0:min(n_max,n_object_found)]:

        n_sample_read = n_sample_read + 1

        #        if (sample.field('Z') > z_min) & (sample.field('Z') < z_max) :
        if 2>1:
            
            x = numpy.double(sample.field('RA'))
            y = numpy.double(sample.field('DEC'))
            index = n_sample_read

            tmp_ang = stomp.IndexedAngularCoordinate(x,y,index,stomp.AngularCoordinate.Equatorial)
            if (x > 0) & (y>-20) & (y<90):
                n_sample_loaded = n_sample_loaded + 1
                sample_map.push_back(tmp_ang)
            #             if (my_map.Contains(tmp_ang) and
            #                 tmp_ang.Eta() > data_eta_min and
            #                 tmp_ang.Eta() < data_eta_max ):
            #                 n_sample_loaded = n_sample_loaded + 1
            #                 sample_map.push_back(tmp_ang)
            # Here I keep all the objects. No masking at all.
                
                
    print colored("%i objects loaded on the map." %n_sample_loaded, 'cyan')
    t1 = time.time()
    # print 'Time spent: ', t1-t0
            
    return sample_in, sample_map


