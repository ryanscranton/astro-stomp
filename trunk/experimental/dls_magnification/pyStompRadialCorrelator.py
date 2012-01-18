from __future__ import division
import numpy
import stomp
import cosmology

"""
A simple, no-frills radial correlation code for python. 
Given one(two) catalogs and a catalog of random points, the code computes 
the auto(cross)-corrlation between the catalog(s) and subtracts off the 
random contribution as in the L&S'93 estimator in log bins of degrees. 
The code can be used as a stand aloneor imported by another function.

Columns for the catalogs are assumed to be:
ra dec weight redshift

Columns for the randoms need only be:
ra dec
The code uses randomized weights and redshifts from the input catalogs and
assigns them to the random points.
"""

def spherical_distance(obj, cat):
    dx = (obj[0] - cat[:,0])*numpy.cos(obj[1]*numpy.pi/180.0)
    dy = (obj[1] - cat[:,1])
    return numpy.sqrt(numpy.power(dx,2)+numpy.power(dy,2))

cosmo = cosmology.MultiEpoch(0.0,5.01)

class RadialCorrelation(object):

    def __init__(self, cat_1, cat_2, rand_1, rand_2=None, n_random=1,
                 r_min=0.001, r_max=10.0, bins_per_decade=5,
                 verbose=False):
        self.verbose = verbose
        self.cat_1 = cat_1[cat_1[:,0].argsort()]
        self.cat_2 = cat_2[cat_2[:,0].argsort()]
        if rand_2 == None:
            while (n_random*len(cat_1) > len(rand_1) or 
                   n_random*len(cat_2) > len(rand_1)):
                print "Requesting more randoms than would be safe"
                print "\tReducing to",n_random-1
                n_random -= 1
            self.n_random = n_random
            self.randoms_1 = numpy.random.permutation(rand_1)[
                :n_random*len(cat_1)]
            self.randoms_2 = numpy.random.permutation(rand_1)[
                :n_random*len(cat_2)]
        else:
            while (n_random*len(cat_1) > len(rand_1) or 
                   n_random*len(cat_2) > len(rand_1)):
                print "Requesting more randoms than would be safe"
                print "\tReducing to",n_random-1
                n_random -= 1
            self.n_random = n_random
            self.randoms_1 = numpy.random.permutation(rand_1)[
                :n_random*len(cat_1)]
            self.randoms_2 = numpy.random.permutation(rand_2)[
                :n_random*len(cat_2)]

        self.r_min = r_min
        self.r_max = r_max
        self.bins_per_decade = bins_per_decade

        self.angular_bins = []
        if bins_per_decade >= 1:
            unit_double = numpy.floor(numpy.log10(r_min))*bins_per_decade
            r = numpy.power(10.0, unit_double/bins_per_decade)
            while r < r_max:
                if r >= self.r_min and r < r_max:
                    self.angular_bins.append(RadialBin(
                            r,numpy.power(10.0,
                                          (unit_double+1.0)/
                                          bins_per_decade)))
                unit_double += 1.0
                r = numpy.power(10.0, unit_double/bins_per_decade)
        else:
            self.angular_bins.append(RadialBin(self.r_min,
                                               self.r_max))
            self.angular_bins[0].r = (self.r_max+self.r_min)/2.0
        if n_random == 0:
            for bin in self.angular_bins:
                bin._gal_rand = 0.0
                bin._rand_gal = 1.0
                bin._rand_rand = 1.0

    def write_result(self, file_name):
        file_out = open(file_name,'w')
        for bin in self.angular_bins:
            if self.n_random > 1:
                bin.rescale_gal_rand(1.0*self.n_random)
                bin.rescale_rand_gal(1.0*self.n_random)
                bin.rescale_rand_rand(1.0*self.n_random)
            file_out.writelines(str(bin.r)+' '+
                                str(bin.compute_correlation())+' '+
                                str(bin.compute_error())+' '+
                                str(bin._gal_gal)+' '+
                                str(bin._gal_rand)+' '+
                                str(bin._rand_gal)+' '+
                                str(bin._rand_rand)+'\n')
        file_out.close()

    def reset_pairs(self):
        for bin in self.angular_bins:
            bin.reset_pairs()

    def find_auto_correlation(self):
        if self.verbose:
            print "Runing Auto-Correlation"
            print "\tComputing Gal-Gal"
        galaxy_tree = stomp.TreeMap(256)
        for obj in self.cat_1:
            if numpy.logical_not(galaxy_tree.AddPoint(
                    stomp.WeightedAngularCoordinate(
                        obj[0],obj[1],obj[2],
                        stomp.AngularCoordinate.Equatorial))):
                print "Failed to add object!"
        for idx,obj in enumerate(self.cat_1):
            if self.verbose and idx%int(self.cat_1.shape[0]/5)==0:
                print "\tRunning Object",idx
            obj_coord = stomp.WeightedAngularCoordinate(
                obj[0],obj[1],obj[2],
                stomp.AngularCoordinate.Equatorial)
            d_A = cosmo.angular_diameter_distance(obj[3])*numpy.pi/180.0
            if self.verbose and idx%int(self.cat_1.shape[0]/5)==0:
                print "\tAngular Min Max"
                print "\t",self.r_min/d_A, self.r_max/d_A
            for bin in self.angular_bins:
                bin.add_to_gal_gal(galaxy_tree.FindWeightedPairs(
                        obj_coord,
                        bin.r_min/d_A,
                        bin.r_max/d_A))
        for rand_it in xrange(self.n_random):
            if self.verbose:
                print "Starting Random Iteration",rand_it+1
            random_tree = stomp.TreeMap(256)
            randoms = self.randoms_1[
                rand_it*len(self.cat_1):(rand_it+1)*len(self.cat_1)]
            weights = numpy.random.permutation(self.cat_1[:,2])
            redshifts = numpy.random.permutation(self.cat_1[:,3])
            randoms = numpy.array([randoms[:,0],randoms[:,1],
                                   weights,redshifts]).transpose()
            if self.verbose:
                print "\tComputing Gal-Rand"
            for idx,rand in enumerate(randoms):
                rand_coord = stomp.WeightedAngularCoordinate(
                    rand[0],rand[1],rand[2],
                    stomp.AngularCoordinate.Equatorial)
                if numpy.logical_not(random_tree.AddPoint(rand_coord)):
                    print "Failed to add object!"
                if self.verbose and idx%int(randoms.shape[0]/5)==0:
                    print "\tRunning Object",idx
                d_A = cosmo.angular_diameter_distance(rand[3])*numpy.pi/180.0
                for bin in self.angular_bins:
                    gal_rand = galaxy_tree.FindWeightedPairs(
                        rand_coord,bin.r_min/d_A,bin.r_max/d_A)
                    bin.add_to_gal_rand(gal_rand)
                    bin.add_to_rand_gal(gal_rand)
            if self.verbose:
                print "\tComputing Rand-Rand"
            for idx,rand in enumerate(randoms):
                if self.verbose and idx%int(randoms.shape[0]/5)==0:
                    print "\tRunning Object",idx
                rand_coord = stomp.WeightedAngularCoordinate(
                    rand[0],rand[1],rand[2],
                    stomp.AngularCoordinate.Equatorial)
                d_A = cosmo.angular_diameter_distance(rand[3])*numpy.pi/180.0
                for bin in self.angular_bins:
                    bin.add_to_rand_rand(random_tree.FindWeightedPairs(
                            rand_coord,bin.r_min/d_A,bin.r_max/d_A))

    def find_cross_correlation(self):
        if self.verbose:
            print "Running Cross Correlation"
            print "\tComputing Gal-Gal"
        galaxy_tree = stomp.TreeMap(256)
        for obj in self.cat_2:
            if numpy.logical_not(galaxy_tree.AddPoint(
                    stomp.WeightedAngularCoordinate(
                        obj[0],obj[1],obj[2],
                        stomp.AngularCoordinate.Equatorial))):
                print "Failed to add object!"
        for idx,obj in enumerate(self.cat_1):
            if self.verbose and idx%int(self.cat_1.shape[0]/5)==0:
                print "\tRunning Object",idx
            obj_coord = stomp.WeightedAngularCoordinate(
                obj[0],obj[1],obj[2],
                stomp.AngularCoordinate.Equatorial)
            d_A = cosmo.angular_diameter_distance(obj[3])*numpy.pi/180.0
            for bin in self.angular_bins:
                bin.add_to_gal_gal(galaxy_tree.FindWeightedPairs(
                        obj_coord,
                        bin.r_min/d_A,
                        bin.r_max/d_A))
        for rand_it in xrange(self.n_random):
            if self.verbose:
                print "Starting Random Iteration",rand_it+1
            rand_1 = self.randoms_1[
                rand_it*len(self.cat_1):(rand_it+1)*len(self.cat_1)]
            weights = numpy.random.permutation(self.cat_1[:,2])
            redshifts = numpy.random.permutation(self.cat_1[:,2])
            rand_1 = numpy.array([rand_1[:,0],rand_1[:,1],
                                  weights,redshifts]).transpose()
            rand_2 = self.randoms_2[
                rand_it*len(self.cat_2):(rand_it+1)*len(self.cat_2)]
            weights = numpy.random.permutation(self.cat_2[:,2])
            redshifts = numpy.random.permutation(self.cat_2[:,2])
            rand_2 = numpy.array([rand_2[:,0],rand_2[:,1],
                                  weights,redshifts]).transpose()
            random_tree = stomp.TreeMap(256)
            for obj in rand_2:
                if numpy.logical_not(random_tree.AddPoint(
                        stomp.WeightedAngularCoordinate(
                            obj[0],obj[1],obj[2],
                            stomp.AngularCoordinate.Equatorial))):
                    print "Failed to add object!"
            if self.verbose:
                print "\tComputing Gal-Rand"
            for idx,obj in enumerate(self.cat_1):
                if self.verbose and idx%int(self.cat_1.shape[0]/5)==0:
                    print "\tRunning Object",idx
                obj_coord = stomp.WeightedAngularCoordinate(
                    obj[0],obj[1],obj[2],
                    stomp.AngularCoordinate.Equatorial)
                d_A = cosmo.angular_diameter_distance(obj[3])*numpy.pi/180.0
                for bin in self.angular_bins:
                    bin.add_to_gal_rand(random_tree.FindWeightedPairs(
                            obj_coord,
                            bin.r_min/d_A,
                            bin.r_max/d_A))
            if self.verbose:
                print "\tComputing Rand-Gal"
            for idx,obj in enumerate(rand_1):
                if self.verbose and idx%int(rand_1.shape[0]/5)==0:
                    print "\tRunning Object",idx
                obj_coord = stomp.WeightedAngularCoordinate(
                    obj[0],obj[1],obj[2],
                    stomp.AngularCoordinate.Equatorial)
                d_A = cosmo.angular_diameter_distance(obj[3])*numpy.pi/180.0
                for bin in self.angular_bins:
                    bin.add_to_rand_gal(galaxy_tree.FindWeightedPairs(
                            obj_coord,
                            bin.r_min/d_A,
                            bin.r_max/d_A))
            if self.verbose:
                print "\tComputing Rand-Rand"
            for idx,obj in enumerate(rand_1):
                if self.verbose and idx%int(rand_1.shape[0]/5)==0:
                    print "\tRunning Object",idx
                obj_coord = stomp.WeightedAngularCoordinate(
                    obj[0],obj[1],obj[2],
                    stomp.AngularCoordinate.Equatorial)
                d_A = cosmo.angular_diameter_distance(obj[3])*numpy.pi/180.0
                for bin in self.angular_bins:
                    bin.add_to_rand_rand(random_tree.FindWeightedPairs(
                            obj_coord,
                            bin.r_min/d_A,
                            bin.r_max/d_A))

class RadialBin(object):
    
    def __init__(self, r_min, r_max):
        self.r_min = r_min
        self.r_max = r_max
        self.r = numpy.power(10.0,0.5*(numpy.log10(r_min)+
                                 numpy.log10(r_max)))
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
    print "Running Unit Test"
    random_1 = numpy.random.uniform(0.01,1.0,size=(10000,4))
    random_2 = numpy.random.uniform(0.01,1.0,size=(10000,4))
    random_3 = numpy.random.uniform(0.01,1.0,size=(1000000,4))
    test_correlator = RadialCorrelation(random_1, random_2, random_3,
                                        n_random=5,verbose=True)
    test_correlator.find_auto_correlation()
    test_correlator.write_result('test_auto_correlation.ascii')
    test_correlator.reset_pairs()
    test_correlator.find_cross_correlation()
    test_correlator.write_result('test_cross_correlation.ascii')
