from __future__ import division
import numpy
import stomp

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
            for bin in self.angular_bins:
                bin.add_to_gal_gal(galaxy_tree.FindWeightedPairs(
                        obj_coord,
                        bin.theta_min,
                        bin.theta_max))
        for rand_it in xrange(self.n_random):
            if self.verbose:
                print "Starting Random Iteration",rand_it+1
            random_tree = stomp.TreeMap(256)
            randoms = numpy.random.permutation(self.randoms)[:len(self.cat_1)]
            weights = numpy.random.permutation(self.cat_1[:,2])
            randoms = numpy.array([randoms[:,0],randoms[:,1],
                                   weights]).transpose()
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
                for bin in self.angular_bins:
                    gal_rand = galaxy_tree.FindWeightedPairs(
                        rand_coord,bin.theta_min,bin.theta_max)
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
                for bin in self.angular_bins:
                    bin.add_to_rand_rand(random_tree.FindWeightedPairs(
                            rand_coord,bin.theta_min,bin.theta_max))

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
            for bin in self.angular_bins:
                bin.add_to_gal_gal(galaxy_tree.FindWeightedPairs(
                        obj_coord,
                        bin.theta_min,
                        bin.theta_max))
        for rand_it in xrange(self.n_random):
            if self.verbose:
                print "Starting Random Iteration",rand_it+1
            rand_1 = numpy.random.permutation(self.randoms)[:len(self.cat_1)]
            weights = numpy.random.permutation(self.cat_1[:,2])
            rand_1 = numpy.array([rand_1[:,0],rand_1[:,1],weights]).transpose()
            rand_2 = numpy.random.permutation(self.randoms)[:len(self.cat_2)]
            weights = numpy.random.permutation(self.cat_2[:,2])
            rand_2 = numpy.array([rand_2[:,0],rand_2[:,1],weights]).transpose()
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
                for bin in self.angular_bins:
                    bin.add_to_gal_rand(random_tree.FindWeightedPairs(
                            obj_coord,
                            bin.theta_min,
                            bin.theta_max))
            if self.verbose:
                print "\tComputing Rand-Gal"
            for idx,obj in enumerate(rand_1):
                if self.verbose and idx%int(rand_1.shape[0]/5)==0:
                    print "\tRunning Object",idx
                obj_coord = stomp.WeightedAngularCoordinate(
                    obj[0],obj[1],obj[2],
                    stomp.AngularCoordinate.Equatorial)
                for bin in self.angular_bins:
                    bin.add_to_rand_gal(galaxy_tree.FindWeightedPairs(
                            obj_coord,
                            bin.theta_min,
                            bin.theta_max))
            if self.verbose:
                print "\tComputing Rand-Rand"
            for idx,obj in enumerate(rand_1):
                if self.verbose and idx%int(rand_1.shape[0]/5)==0:
                    print "\tRunning Object",idx
                obj_coord = stomp.WeightedAngularCoordinate(
                    obj[0],obj[1],obj[2],
                    stomp.AngularCoordinate.Equatorial)
                for bin in self.angular_bins:
                    bin.add_to_rand_rand(random_tree.FindWeightedPairs(
                            obj_coord,
                            bin.theta_min,
                            bin.theta_max))

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
    print "Running Unit Test"
    random_1 = numpy.random.uniform(0,1,size=(10000,3))
    random_2 = numpy.random.uniform(0,1,size=(10000,3))
    random_3 = numpy.random.uniform(0,1,size=(1000000,3))
    test_correlator = AngularCorrelation(random_1, random_2, random_3,
                                         n_random=5,verbose=True)
    #test_correlator.find_auto_correlation()
    #test_correlator.write_result('test_auto_correlation.ascii')
    test_correlator.reset_pairs()
    test_correlator.find_cross_correlation()
    test_correlator.write_result('test_cross_correlation.ascii')
