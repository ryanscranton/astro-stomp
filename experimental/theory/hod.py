import numpy
import param
from scipy import interpolate

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

    def first_moment(self, mass, z):
        return 1.0

    def second_moment(self, mass, z):
        return 1.0

    def nth_moment(self, mass, z, n=3):
        exp_nth = 0.0
        if n in self._hod:
            exp_nth = self._hod[n](mass, z)
        else:
            exp_nth = self.first_moment(mass, z)**n
            alpha_m2 = self.second_moment(mass, z)/self.first_moment(mass, z)**2
            for j in xrange(n):
                exp_nth *= (j*self.alpha_m_2(mass, z) - j + 1)
        return exp_nth

class HODPoisson(HaloModel):
    def __init__(self):
        pass

class HODBinomial(HaloModel):
    def __init__(self, n_max, min_mass, mass_max, p_m_spline):
        pass

class HODKravtsov(HaloModel):
    def __init__(self, halo_param=None, **kws):
        HaloModel.__init__()

        if halo_param=None:
            halo_param = param.HaloModelParams(**kws)
        self.min_mass = halo_param.m11_msunh
        self.mass_one = halo_param.m13_msunh
        self.mass_cut = halo_param.mcut_msunh
        self.alpha = halo_param.alphaexp

        self._hod[1] = self.first_moment
        self._hod[2] = self.second_moment

    def first_moment(self, mass, z):
        return (self.central_first_moment(mass, z) +
                self.satellite_first_moment(mass, z))

    def second_moment(self, mass, z):
        return self.satellite_second_moment(mass, z)

    def central_first_moment(self, mass, z):
        if mass >= self.min_mass:
            return 1.0
        else:
            return 0.0

    def satellite_first_moment(self, mass, z):
        return mass/self.mass_one*numpy.exp(-self.mass_cut/mass)

    def second_moment(self, mass, z):
        return self.satellite_second_moment(mass, z)

    def satellite_second_moment(self, mass, z):
        return self.satellite_first_moment(mass, z)**2
