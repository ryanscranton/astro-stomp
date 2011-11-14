from __future__ import division
import numpy

"""
A simple, no-frills correlation code for python. Given one(two) catalogs and a
catalog of random points, the code computes the auto(cross)-corrlation between
the catalog(s) and subtracts off the random contribution as in the L&S'93
estimator in log bins of degrees. The code can be used as a stand alone
or imported by another function.
"""

def spherical_distance(obj, cat):
    dx = (obj[0] - cat[:,0])*numpy.cos(obj[1]*numpy.pi/180.0)
    dy = (obj[1] - cat[:,1])
    return numpy.sqrt(numpy.power(dx,2)+numpy.power(dy,2))

class AngularCorrelation(object):

    def __init__(self, cat_1, cat_2, rand, n_random=1,
                 theta_min=0.001, theta_max=0.5, bins_per_decade=5,
                 verbose=False):
        self.verbose = verbose
        while (n_random*len(cat_1) > len(rand) or
               n_random*len(cat_2) > len(rand)):
            print "Requesting more randoms than would be safe"
            print "\tReducing to",n_random-1
            n_random -= 1
        self.n_random = n_random
        self.cat_1 = cat_1[cat_1[:,0].argsort()]
        self.cat_2 = cat_2[cat_2[:,0].argsort()]
        self.randoms = rand

        self.theta_min = theta_min
        self.theta_max = theta_max
        self.bins_per_decade = bins_per_decade

        self.angular_bins = []
        if bins_per_decade >= 1:
            unit_double = numpy.floor(numpy.log10(theta_min))*bins_per_decade
            theta = numpy.power(10.0, unit_double/bins_per_decade)
            while theta < theta_max:
                if theta >= self.theta_min and theta < theta_max:
                    self.angular_bins.append(AngularBin(
                            theta,numpy.power(10.0,
                                              (unit_double+1.0)/
                                              bins_per_decade)))
                unit_double += 1.0
                theta = numpy.power(10.0, unit_double/bins_per_decade)
        else:
            self.angular_bins.append(AngularBin(self.theta_min,
                                                self.theta_max))
            self.angular_bins[0].theta = (self.theta_max+self.theta_min)/2.0

    def write_result(self, file_name):
        file_out = open(file_name,'w')
        for bin in self.angular_bins:
            if self.n_random > 1:
                bin.rescale_gal_rand(1.0*self.n_random)
                bin.rescale_rand_gal(1.0*self.n_random)
                bin.rescale_rand_rand(1.0*self.n_random)
            file_out.writelines(str(bin.theta)+' '+
                                str(bin.compute_correlation())+' '+
                                str(bin.compute_error())+' '+
                                str(bin._gal_gal)+' '+
                                str(bin._gal_rand)+' '+
                                str(bin._rand_gal)+' '+
                                str(bin._rand_rand)+'\n')
        file_out.close()

    def find_auto_correlation(self):
        if self.verbose:
            print "Runing Auto-Correlation"
            print "\tComputing Gal-Gal"
        for idx,obj in enumerate(self.cat_1):
            if self.verbose and idx%int(self.cat_1.shape[0]/5)==0:
                print "\tRunning Object",idx
            ra_max = self.theta_max/numpy.cos(obj[1]*numpy.pi/180.0)
            tmp_cat = self.cat_1[numpy.logical_and(
                    numpy.logical_and(
                        self.cat_1[:,0]>obj[0]-ra_max,
                        self.cat_1[:,0]<obj[0]+ra_max),
                    numpy.logical_and(
                        self.cat_1[:,1]>obj[1]-self.theta_max,
                        self.cat_1[:,1]<obj[1]+self.theta_max))]
            dist = spherical_distance(obj, tmp_cat)
            for bin in self.angular_bins:
                gal_gal_mask = numpy.logical_and(dist>bin.theta_min, 
                                                 dist<bin.theta_max)
                bin.add_to_gal_gal(numpy.sum(obj[2]*
                                             tmp_cat[gal_gal_mask,2]))
        for rand_it in xrange(self.n_random):
            if self.verbose:
                print "Starting Random Itteration",rand_it+1
            rand = numpy.random.permutation(self.randoms)[:len(self.cat_1)]
            weights = numpy.random.permutation(self.cat_1[:,2])
            rand = numpy.array([rand[:,0],rand[:,1],weights]).transpose()
            if self.verbose:
                print "\tComputing Gal-Rand"
            for idx,obj in enumerate(self.cat_1):
                if self.verbose and idx%int(self.cat_1.shape[0]/5)==0:
                    print "\tRunning Object",idx
                ra_max = self.theta_max/numpy.cos(obj[1]*numpy.pi/180.0)
                tmp_rand = rand[numpy.logical_and(
                        numpy.logical_and(
                            rand[:,0]>obj[0]-ra_max,
                            rand[:,0]<obj[0]+ra_max),
                        numpy.logical_and(
                            rand[:,1]>obj[1]-self.theta_max,
                            rand[:,1]<obj[1]+self.theta_max))]
                dist = spherical_distance(obj, tmp_rand)
                for bin in self.angular_bins:
                    gal_rand_mask = numpy.logical_and(dist>bin.theta_min, 
                                                      dist<bin.theta_max)
                    gal_rand = numpy.sum(
                        obj[2]*tmp_rand[gal_rand_mask,2])
                    bin.add_to_gal_rand(gal_rand)
                    bin.add_to_rand_gal(gal_rand)
            if self.verbose:
                print "\tComputing Rand-Rand"
            for idx,obj in enumerate(rand):
                if self.verbose and idx%int(rand.shape[0]/5)==0:
                    print "\tRunning Object",idx
                ra_max = self.theta_max/numpy.cos(obj[1]*numpy.pi/180.0)
                tmp_rand = rand[numpy.logical_and(
                        numpy.logical_and(
                            rand[:,0]>obj[0]-ra_max,
                            rand[:,0]<obj[0]+ra_max),
                        numpy.logical_and(
                            rand[:,1]>obj[1]-self.theta_max,
                            rand[:,1]<obj[1]+self.theta_max))]
                dist = spherical_distance(obj, tmp_rand)
                for bin in self.angular_bins:
                    rand_rand_mask = numpy.logical_and(dist>bin.theta_min, 
                                                       dist<bin.theta_max)
                    bin.add_to_rand_rand(numpy.sum(
                            obj[2]*tmp_rand[rand_rand_mask,2]))

    def find_cross_correlation(self):
        if self.verbose:
            print "Running Cross Correlation"
            print "\tComputing Gal-Gal"
        for idx,obj in enumerate(self.cat_1):
            if self.verbose and idx%int(self.cat_1.shape[0]/5)==0:
                print "\tRunning Object",idx
            ra_max = self.theta_max/numpy.cos(obj[1]*numpy.pi/180.0)
            tmp_cat = self.cat_2[numpy.logical_and(
                    numpy.logical_and(
                        self.cat_2[:,0]>obj[0]-ra_max,
                        self.cat_2[:,0]<obj[0]+ra_max),
                    numpy.logical_and(
                        self.cat_2[:,1]>obj[1]-self.theta_max,
                        self.cat_2[:,1]<obj[1]+self.theta_max))]
            dist = spherical_distance(obj, tmp_cat)
            for bin in self.angular_bins:
                gal_gal_mask = numpy.logical_and(dist>bin.theta_min, 
                                                 dist<bin.theta_max)
                bin.add_to_gal_gal(numpy.sum(obj[2]*
                                             tmp_cat[gal_gal_mask,2]))
        for rand_it in xrange(self.n_random):
            if self.verbose:
                print "Starting Random Itteration",rand_it+1
            rand_1 = numpy.random.permutation(self.randoms)[:len(self.cat_1)]
            weights = numpy.random.permutation(self.cat_1[:,2])
            rand_1 = numpy.array([rand_1[:,0],rand_1[:,1],weights]).transpose()
            rand_2 = numpy.random.permutation(self.randoms)[:len(self.cat_2)]
            weights = numpy.random.permutation(self.cat_2[:,2])
            rand_2 = numpy.array([rand_2[:,0],rand_2[:,1],weights]).transpose()
            if self.verbose:
                print "\tComputing Gal-Rand"
            for idx,obj in enumerate(self.cat_1):
                if self.verbose and idx%int(self.cat_1.shape[0]/5)==0:
                    print "\tRunning Object",idx
                ra_max = self.theta_max/numpy.cos(obj[1]*numpy.pi/180.0)
                tmp_rand = rand_2[numpy.logical_and(
                        numpy.logical_and(
                            rand_2[:,0]>obj[0]-ra_max,
                            rand_2[:,0]<obj[0]+ra_max),
                        numpy.logical_and(
                            rand_2[:,1]>obj[1]-self.theta_max,
                            rand_2[:,1]<obj[1]+self.theta_max))]
                dist = spherical_distance(obj, tmp_rand)
                for bin in self.angular_bins:
                    gal_rand_mask = numpy.logical_and(dist>bin.theta_min, 
                                                      dist<bin.theta_max)
                    bin.add_to_gal_rand(numpy.sum(
                            obj[2]*tmp_rand[gal_rand_mask,2]))
            if self.verbose:
                print "\tComputing Rand-Gal"
            for idx,obj in enumerate(rand_1):
                if self.verbose and idx%int(rand_1.shape[0]/5)==0:
                    print "\tRunning Object",idx
                ra_max = self.theta_max/numpy.cos(obj[1]*numpy.pi/180.0)
                tmp_cat = self.cat_2[numpy.logical_and(
                        numpy.logical_and(
                            self.cat_2[:,0]>obj[0]-ra_max,
                            self.cat_2[:,0]<obj[0]+ra_max),
                        numpy.logical_and(
                            self.cat_2[:,1]>obj[1]-self.theta_max,
                            self.cat_2[:,1]<obj[1]+self.theta_max))]
                dist = spherical_distance(obj, tmp_cat)
                for bin in self.angular_bins:
                    rand_gal_mask = numpy.logical_and(dist>bin.theta_min, 
                                                      dist<bin.theta_max)
                    bin.add_to_rand_gal(numpy.sum(
                            obj[2]*tmp_cat[rand_gal_mask,2]))
            if self.verbose:
                print "\tComputing Rand-Rand"
            for idx,obj in enumerate(rand_1):
                if self.verbose and idx%int(rand_1.shape[0]/5)==0:
                    print "\tRunning Object",idx
                ra_max = self.theta_max/numpy.cos(obj[1]*numpy.pi/180.0)
                tmp_rand = rand_2[numpy.logical_and(
                        numpy.logical_and(
                            rand_2[:,0]>obj[0]-ra_max,
                            rand_2[:,0]<obj[0]+ra_max),
                        numpy.logical_and(
                            rand_2[:,1]>obj[1]-self.theta_max,
                            rand_2[:,1]<obj[1]+self.theta_max))]
                dist = spherical_distance(obj, tmp_rand)
                for bin in self.angular_bins:
                    rand_rand_mask = numpy.logical_and(dist>bin.theta_min, 
                                                       dist<bin.theta_max)
                    bin.add_to_rand_rand(numpy.sum(
                            obj[2]*tmp_rand[rand_rand_mask,2]))

class AngularBin(object):
    
    def __init__(self, theta_min, theta_max):
        self.theta_min = theta_min
        self.theta_max = theta_max
        self.theta = numpy.power(10.0,0.5*(numpy.log10(theta_min)+
                                 numpy.log10(theta_max)))
        self.reset_pairs()

    def reset_pairs(self):
        self._gal_gal = 0.0
        self._gal_rand = 0.0
        self._rand_gal = 0.0
        self._rand_rand = 0.0

    def add_to_gal_gal(self, pairs):
        self._gal_gal += pairs

    def add_to_rand_gal(self, pairs):
        self._rand_gal += pairs

    def add_to_gal_rand(self, pairs):
        self._gal_rand += pairs

    def add_to_rand_rand(self, pairs):
        self._rand_rand += pairs

    def rescale_gal_rand(self, weight):
        self._gal_rand = self._gal_rand/weight

    def rescale_rand_gal(self, weight):
        self._rand_gal = self._rand_gal/weight

    def rescale_rand_rand(self, weight):
        self._rand_rand = self._rand_rand/weight

    def compute_correlation(self):
        return (self._gal_gal-self._gal_rand-
                self._rand_gal+self._rand_rand)/self._rand_rand

    def compute_error(self):
        return 1.0/numpy.sqrt(self._gal_gal)

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f","--fore",dest="fore",default="",
                      action="store",type="str",
                      help="Foreground file to load")
    parser.add_option("-b","--back",dest="back",default="",
                      action="store",type="str",help="Background File to load")
    parser.add_option("-r","--rand",dest="rand",default="",
                      action="store",type="str",help="Randoms File To Load")
    parser.add_option("-o","--output_file",dest="output_file",
                      default="Wtheta_test",action="store",type="str",
                      help="Output file name")
    parser.add_option("-a","--auto_corr",dest="auto_corr",
                      default="false",action="store",type="str",
                      help="Run Auto-Correlation")
    parser.add_option("-n","--n_randoms",dest="n_randoms",default=1,
                      action="store",type="int",help="# Randoms")
    parser.add_option('--theta_min',dest='theta_min',default=0.001,
                      action="store",type="float",help="Minimum Angle")
    parser.add_option('--theta_max',dest='theta_max',default=0.5,
                      action="store",type="float",help="Maximum Angle")
    parser.add_option('--bins_per_decade',dest='bins_per_decade',default=5,
                      action="store",type="int",help="# Bins Per Decade")
    (options, args) = parser.parse_args()

    fore = numpy.loadtxt(options.fore)
    back = numpy.loadtxt(options.back)
    rand = numpy.loadtxt(options.rand)
    rand = rand[numpy.logical_and(
            numpy.logical_and(rand[:,2]>1000,rand[:,2]<9000),
            numpy.logical_and(rand[:,3]>1000,rand[:,3]<9000))]

    wtheta = AngularCorrelation(fore, back, rand, options.n_randoms,
                                theta_min=options.theta_min,
                                theta_max=options.theta_max, 
                                bins_per_decade=options.bins_per_decade)
    if options.auto_corr == "false":
        wtheta.find_cross_correlation()
    else:
        wtheta.find_auto_correlation()
    wtheta.write_result(options.output_file)
