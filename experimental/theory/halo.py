from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import integrate
import param  # cosmological parameter object from Cosmopy
import numpy
import cosmology
import mass_function
import camb
import hod

"""Classes for describing a basic halo model.

Given an input HOD object and cosmological and halo parameters, a halo object
should be able to generate a non-linear power spectrum for dark matter,
galaxies or their cross-spectrum.  It should also be able to return the halo
profile as a function of mass and radius and its Fourier transform as a
function of mass and wavenumber.
"""

__author__ = "Ryan Scranton <ryan.scranton@gmail.com"

class Halo(object):
    """Basic halo model object.

    Given an HOD object and cosmological and halo parameters, a halo object
    should be able to generate a non-linear power spectrum for dark matter,
    galaxies or their cross-spectrum.  It should also be able to return
    the halo profile as a function of mass and radius and its Fourier transform
    as a function of mass and wavenumber.
    """
    def __init__(self, hod, redshift=None, camb_param=None,
                 halo_param=None, **kws):
        if camb_param is None:
            camb_param = param.CambParams(**kws)

        if halo_param is None:
            halo_param = param.HaloModelParams(**kws)

        if redshift is None: redshift = 0.0
        self.redshift = redshift
        z_min = 0.0
        z_max = self.redshift + 1.0

        self.cosmo = cosmology.Cosmology(z_min, z_max, camb_param)

        self.mass = mass_function.MassFunction(
            self.redshift, camb_param, halo_param)

        self.hod = hod

        self.camb = camb.CambWrapper(camb_param)
        self.camb.set_redshift(redshift)
        self.camb.run()

        # self.nbar
        self.calculate_nbar()

    def linear_power(self, k):
        return self.camb.linear_power(k)

    def power_mm(self, k):
        return self._Phh_mm(k) + self._PP_mm(k)

    def power_gm(self, k):
        return self._Phh_gm(k) + self._PP_gm(k)

    def power_gg(self, k):
        return self._Phh_gg(k) + self._PP_gg(k)

    def _Phh_mm(self, k):
        pass

    def _PP_mm(self, k):
        pass

    def _Phh_gm(self, k):
        pass

    def _PP_gm(self, k):
        pass

    def _Phh_gg(self, k):
        pass

    def _PP_gg(self, k):
        pass

    def y(self, k, mass):
        pass

    
