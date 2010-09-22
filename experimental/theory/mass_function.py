from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import integrate
import param  # cosmological parameter object from Cosmopy
import numpy
import cosmology

"""Classes for encoding basic cosmological parameters and quantities.

Given a cosmology, we can derive a mass function object.  This should be able
to return the abundance of halo at a given mass as well as the mean halo bias.
Likewise, these quantities should be available as a function of nu:

nu = (delta_c/sigma(M))^2

This is the ratio of the over-density necessary for linear, spherical
collapse and the RMS density fluctuations as a function of mass.  M*, the mass
at which nu == 1, represents the break between halos that are mostly collapsed
and those where fluctuations large enough for collapse are exponentially
suppressed.
"""

__author__ = "Ryan Scranton <ryan.scranton@gmail.com"

class MassFunction(object):
    """Object representing a mass function for a given input cosmology.

    A MassFunction object can return a properly normalized halo abundance or
    halo bias as a function of halo mass or as a function of nu, as well as
    translate between mass and nu.
    """
    def __init__(self, redshift=0.0, camb_param=None, halo_param=None, **kws):
        # Hard coded, but we shouldn't expect halos outside of this range.
        mass_min = 1.0e9
        mass_max = 5.0e17

        dlog_mass = (numpy.log10(mass_max) - numpy.log10(mass_min))/100
        self.log_mass_max = numpy.log10(mass_max) + dlog_mass
        self.log_mass_min = numpy.log10(mass_min) - dlog_mass

        self._log_mass_array = numpy.arange(
            self.log_mass_min, self.log_mass_max + dlog_mass, dlog_mass)

        if camb_param is None:
            camb_param = param.CambParams(**kws)

        if halo_param is None:
            halo_param = param.HaloModelParams(**kws)

        self.stq = halo_param.stq
        self.st_little_a = halo_param.st_little_a

        self.redshift = redshift
        z_min = 0.0
        z_max = self.redshift + 1.0

        self.cosmo = cosmology.Cosmology(z_min, z_max, camb_param)

        self._initialize_splines()
        self._normalize()

    def _initialize_splines(self):
        self._nu_array = numpy.zeros_like(self._log_mass_array)

        for idx in xrange(self._log_mass_array.size):
            mass = 10.0**self._log_mass_array[idx]
            self._nu_array[idx] = self.cosmo.nu_m(mass, self.redshift)

        self.nu_min = 1.001*self._nu_array[0]
        self.nu_max = 0.999*self._nu_array[-1]

        self._nu_spline = InterpolatedUnivariateSpline(
            self._log_mass_array, self._nu_array)
        self._log_mass_spline = InterpolatedUnivariateSpline(
            self._nu_array, self._log_mass_array)

    def _normalize(self):
        self.f_norm = 1.0
        norm, norm_err = integrate.quad(self.f_nu, self.nu_min, self.nu_max)
        self.f_norm = 1.0/norm

        self.bias_norm = 1.0
        norm, norm_err = integrate.quad(lambda x: self.f_nu(x)*self.bias_nu(x),
                                        self.nu_min, self.nu_max)
        self.bias_norm = 1.0/norm

    def f_nu(self, nu):
        nu_prime = nu*self.st_little_a
        return (
            self.f_norm*(1.0 + nu_prime**(-1.0*self.stq))*
            numpy.sqrt(nu_prime)*numpy.exp(-0.5*nu_prime)/nu)

    def f_m(self, mass):
        return self.f_nu(self._nu_spline(numpy.log10(mass)))

    def bias_nu(self, nu, redshift=None):
        nu_prime = nu*self.st_little_a
        delta_c = self.cosmo.delta_c(redshift)
        return self.bias_norm*(
            1.0 + (nu - 1.0)/delta_c +
            2.0*self.stq/(delta_c*(1.0 + nu_prime**self.stq)))

    def bias_m(self, mass, redshift=None):
        return self.bias_nu(self._nu_spline(numpy.log10(mass)))

    def nu(self, mass):
        return self._nu_spline(numpy.log10(mass))

    def mass(self, nu):
        return 10.0**self._log_mass_spline(nu)

    def m_star(self):
        return 10.0**self._log_mass_spline(1.0)
