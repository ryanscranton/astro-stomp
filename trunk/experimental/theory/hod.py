import numpy
import param
from scipy import interpolate
from scipy import special

"""Classes for implementing a halo occupation distribution.

A halo occupation distribution describes the population of a given dark matter
halo.  In principle, this distribution can be a function of any number of halo
and galaxy variables.  For the purposes of this implementation, we allow for
the HOD to vary as a function of halo mass and redshift.  Any HOD-derived
object should be capable of returning the first and second moments of the
distribution as a function of mass and redshift.
"""

__author__ = "Chris Morrison <morrison.chrisb@gmail.com>"

class HOD(object):
    """Base class for a halo occupation distribution.

    The exact form of an HOD will vary depending on the galaxy type and
    redshift.  Our base class defines the API that all of the derived class
    instances need to implement.
    """
    def __init__(self):
        self._hod = {}
        self._hod[1] = self.first_moment
        self._hod[2] = self.second_moment

    def first_moment(self, mass, z=None):
        return 1.0

    def second_moment(self, mass, z=None):
        return 1.0

    def nth_moment(self, mass, z=None, n=3):
        exp_nth = 0.0
        if n in self._hod:
            exp_nth = self._hod[n](mass, z)
        else:
            exp_nth = self.first_moment(mass, z)**n
            alpha_m2 = self.second_moment(mass, z)/self.first_moment(mass, z)**2
            for j in xrange(n):
                exp_nth *= (j*alpha_m2 - j + 1)
        return exp_nth

    def set_halo(halo_param):
        pass

    def write(self, output_file_name):
        mass_max = 1.0e16
        mass_min = 1.0e9

        dln_mass = (numpy.log(mass_max) - numpy.log(mass_min))/200
        ln_mass_max = numpy.log(mass_max) + dln_mass
        ln_mass_min = numpy.log(mass_min) - dln_mass

        ln_mass_array = numpy.arange(
            ln_mass_min, ln_mass_max + dln_mass, dln_mass)

        f = open(output_file_name, "w")
        for ln_mass in ln_mass_array:
            mass = numpy.exp(ln_mass)
            f.write("%1.10f %1.10f %1.10f %1.10f\n" % (
                mass, self.first_moment(mass), self.second_moment(mass),
                self.nth_moment(mass,None,3)))
        f.close()

class HODPoisson(HOD):

    def __init__(self):
        pass

class HODBinomial(HOD):

    def __init__(self, n_max, min_mass, mass_max, p_m_spline):
        pass

class HODKravtsov(HOD):
    def __init__(self, halo_param=None, **kws):
        HOD.__init__(self)

        if halo_param is None:
            halo_param = param.HaloModelParams(**kws)
        self.set_halo(halo_param)

    def set_halo(self, halo_param):
        self.min_mass = halo_param.m11_msunh
        self.mass_one = halo_param.m13_msunh
        self.mass_cut = halo_param.mcut_msunh
        self.alpha = halo_param.alphaexp

    def first_moment(self, mass, z=None):
        return (self.central_first_moment(mass) +
                self.satellite_first_moment(mass))

    def second_moment(self, mass, z=None):
        return self.satellite_second_moment(mass)

    def central_first_moment(self, mass):
        if mass >= self.min_mass:
            return 1.0
        else:
            return 0.0

    def satellite_first_moment(self, mass):
        return mass/self.mass_one*numpy.exp(-self.mass_cut/mass)

    def satellite_second_moment(self, mass):
        return (self.satellite_first_moment(mass))**2

class HODZheng(HOD):

    def __init__(self, M_min=10**11.83, sigma=0.30, 
                 M_0=10**11.53, M_1=10**13.02, alpha=1.0, w=1.0):
        self.M_min = M_min
        self.sigma = sigma
        self.M_0 = M_0
        self.M_1 = M_1
        self.alpha = alpha
        self.w = w
        HOD.__init__(self)

    def first_moment(self, mass, z=None):
        return self.central_term(mass)*(
            1+self.satellite_term(mass))

    def second_moment(self, mass, z=None):
        return self.w*(self.central_term(mass)*self.satellite_term(mass))**2

    def central_term(self, mass):
        return 0.5*(1+special.erf((numpy.log10(mass) - numpy.log10(self.M_min))/
                                  self.sigma))

    def satellite_term(self, mass):
        diff = mass - self.M_0
        if diff >=0.0:
            return numpy.power((mass - self.M_0)/self.M_1,self.alpha)
        else:
            return 0.0

class HODMandelbaum(HOD):

    def __init__(self, M0=1.0e13, norm=1.0e13):
        HOD.__init__(self)

        self.M0 = M0
        self.norm = norm

        self._hod[1] = self.first_moment
        self._hod[2] = self.second_moment

    def first_moment(self, mass, z=None):
        if mass > self.M0:
            return mass/self.norm
        if mass <= self.M0:
            return mass*mass/(self.M0*self.norm)

    def second_moment(self, mass, z=None):
        return self.first_moment(mass)**2
