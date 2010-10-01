import numpy
import param
from scipy import interpolate

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

        self._hod[1] = self.first_moment
        self._hod[2] = self.second_moment

    def set_halo(halo_param):
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
        return self.satellite_first_moment(mass)**2
