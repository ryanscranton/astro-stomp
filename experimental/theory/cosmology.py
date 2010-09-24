from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import integrate
import param  # cosmological parameter object from Cosmopy
import numpy

"""Classes for encoding basic cosmological parameters and quantities.

This is inspired by the CosmoPy package.  The intent is to have a singleton
object that contains all of the basic cosmological parameters as well as those
expected by the CAMB software.  In addition, we want a single class object
that can take those parameters and return standard cosmological quantities
(comoving distance, linear growth factor, etc.) as a function of redshift.
"""

__author__ = "Ryan Scranton <ryan.scranton@gmail.com>"


class Cosmology(object):
    """Container class for calculating cosmological values at a given redshift.

    Given an initial set of redshift bounds and a set of parameters, Cosmology
    can return the comoving distance, angular diameter distance or luminosity
    distance for any redshift over the input interval.  Likewise, it can return
    a redshift, given an input comoving distance, as well as the growth factor
    D(z).
    """
    def __init__(self, z_min, z_max, camb_param=None, **kws):
        dz = (z_max - z_min)/100
        self.z_max = z_max + dz
        self.z_min = z_min - dz
        if self.z_min < 0.0: self.z_min = 0.0


        dz = (self.z_max - self.z_min)/100
        self._z_array = numpy.arange(self.z_min, self.z_max + dz, dz)
        self._chi_array = numpy.zeros_like(self._z_array)
        self._growth_array = numpy.zeros_like(self._z_array)

        if camb_param is None:
            camb_param = param.CambParams(**kws)

        self.omega_m0 = camb_param.omega_cdm + camb_param.omega_baryon
        self.omega_b0 = camb_param.omega_baryon
        self.omega_l0 = camb_param.omega_lambda
        self.h = camb_param.hubble/100.0
        self.H0 = self.h/3000.0
        self.omega_r0 = 4.17e-5/(self.h**2)

        self.k_min = 0.001
        self.k_max = 100.0
        self.n = camb_param.scalar_spectral_index[0]
        self.delta_H = (
            1.94e-5*self.omega_m0**(-0.785 - 0.05*numpy.log(self.omega_m0))*
            numpy.exp(-0.95*(self.n - 1) - 0.169*(self.n - 1)**2))

        self._initialize_splines()

    def _initialize_splines(self):
        for idx in xrange(self._z_array.size):
            dist, dist_err = integrate.quad(self.E, 0.0, self._z_array[idx])
            self._chi_array[idx] = dist
        self._chi_spline = InterpolatedUnivariateSpline(self._z_array,
                                                        self._chi_array)
        self._z_spline = InterpolatedUnivariateSpline(self._chi_array,
                                                      self._z_array)

        a_min = 1.0e-5
        self.growth_norm, growth_norm_err = integrate.quad(
            self._growth_integrand, a_min, 1.0)
        self.growth_norm *= 2.5*self.omega_m0/(self.E(0.0)*self.H0)

        for idx in xrange(self._z_array.size):
            a = 1.0/(1.0 + self._z_array[idx])
            growth, growth_err = integrate.quad(
                self._growth_integrand, a_min, a)
            growth *= 2.5*self.omega_m0/(self.E(self._z_array[idx])*self.H0)
            self._growth_array[idx] = growth/self.growth_norm

        self._growth_spline = InterpolatedUnivariateSpline(self._z_array,
                                                           self._growth_array)

    def E(self, redshift):
        return 1.0/(self.H0*numpy.sqrt(
            self.omega_l0 + self.omega_m0*(1 + redshift)**3 +
            (1 - self.omega_l0 - self.omega_m0)*(1 + redshift)**2))

    def _growth_integrand(self, redshift):
        return ((1.0 + redshift)*self.E(redshift)*self.H0)**3

    def flat(self):
        omega_total = self.omega_m0 + self.omega_l0
        if omega_total <= 1.0001 and omega_total >= 0.9999:
            return True
        else:
            return False

    def open(self):
        omega_total = self.omega_m0 + self.omega_l0
        if omega_total <= 0.9999:
            return True
        else:
            return False

    def closed(self):
        omega_total = self.omega_m0 + self.omega_l0
        if omega_total > 1.0001:
            return True
        else:
            return False

    def _reset_redshift_maximum(self, redshift):
        dz = (redshift - self.z_min)/100
        self.z_max = redshift
        self._z_array = numpy.arange(self.z_min, self.z_max + dz, dz)
        self._chi_array = numpy.zeros_like(self._z_array)
        self._growth_array = numpy.zeros_like(self._z_array)

        self._initialize_splines()

    def _reset_redshift_minimum(self, redshift):
        self.z_min = redshift
        dz = (self.z_max - self.z_min)/100
        self._z_array = numpy.arange(self.z_min, self.z_max+dz, dz)
        self._chi_array = numpy.zeros_like(self._z_array)
        self._growth_array = numpy.zeros_like(self._z_array)

        self._initialize_splines()

    def comoving_distance(self, redshift):
        distance = 0.0

        if redshift < self.z_max and redshift > self.z_min:
            distance = self._chi_spline(redshift)[0]
        else:
            if redshift > self.z_max:
                self._reset_redshift_maximum(redshift)
                distance = self.comoving_distance(redshift)
            else:
                if redshift > 0.00001:
                    self._reset_redshift_minimum(redshift)
                    distance = self.comoving_distance(redshift)

        return distance

    def luminosity_distance(self, redshift):
        return (1.0 + redshift)*self.comoving_distance(redshift)

    def angular_diameter_distance(self, redshift):
        return self.comoving_distance(redshift)/(1.0 + redshift)

    def redshift(self, comoving_distance):
        return self._z_spline(comoving_distance)[0]

    def growth_factor(self, redshift):
        """Linear growth factor, normalized to unity at z = 0."""
        growth = 1.0

        if redshift < self.z_max and redshift > self.z_min:
            growth = self._growth_spline(redshift)[0]
        else:
            if redshift > self.z_max:
                self._reset_redshift_maximum(redshift)
                growth = self.growth_factor(redshift)
            else:
                if redshift > 0.00001:
                    self._reset_redshift_minimum(redshift)
                    growth = self.growth_factor(redshift)

        return growth

    def omega_m(self, redshift):
        if redshift is None: redshift = 0.0
        return self.omega_m0*(1.0 + redshift)**3/(
            self.omega_l0 + self.omega_r0*(1.0 + redshift)**2 +
            self.omega_m0*(1.0 + redshift)**3)

    def omega_l(self, redshift):
        if redshift is None: redshift = 0.0
        return self.omega_l0/(
            self.omega_l0 + self.omega_r0*(1.0 + redshift)**2 +
            self.omega_m0*(1.0 + redshift)**3)

    def delta_c(self, redshift=None):
        """Over-density threshold for linear, spherical collapse."""
        # Fitting function taken from NFW97
        delta_c = 0.15*(12.0*numpy.pi)**(2.0/3.0)
        if self.open():
            delta_c *= self.omega_m(redshift)**0.0185
        if self.flat() and self.omega_m0 < 1.0001:
            delta_c *= self.omega_m(redshift)**0.0055

        return delta_c/self.growth_factor(redshift)

    def delta_v(self, redshift=None):
        """Over-density for a collapsed, virialized halo."""
        # Fitting function taken from NFW97
        delta_v = 178.0
        if self.open():
            delta_v /= self.omega_m(redshift)**0.7
        if self.flat() and self.omega_m0 < 1.0001:
            delta_v /= self.omega_m(redshift)**0.55

        return delta_v/self.growth_factor(redshift)

    def rho_crit(self, redshift):
        """Critical density in solar masses per cubic Mpc."""
        rho_crit = (self.omega_l0 + self.omega_r0*(1.0 + redshift)**2 +
                    self.omega_m0*(1.0 + redshift)**3)
        return rho_crit*1.0e-29*1.0e-33*2.937999e+73

    def rho_bar(self, redshift):
        """Matter density in solar masses per cubic Mpc."""
        return self.rho_crit(redshift)*self.omega_m(redshift)

    def _bbks(self, k):
        """BBKS transfer function."""
        q = k/self.omega_m0/self.h**2*numpy.exp(
            self.omega_b0 + numpy.sqrt(2*self.h)*self.omega_b0/self.omega_m0)
        return (numpy.log(1.0 + 2.34*q)/(2.34*q)*
                (1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**(-0.25))

    def delta_k(self, k, redshift=None):
        """k^2*P(k)/2*pi^2: dimensionless linear power spectrum.

        If redshift is omitted the results are given for z = 0.
        """
        delta_k = (
            self.delta_H**2*(3000.0*k/self.h)**(3 + self.n)*self._bbks(k)**2)
        if not redshift is None:
            delta_k *= self.growth_factor(redshift)**2
        return delta_k

    def linear_power(self, k, redshift=None):
        return 2.0*numpy.pi**2*self.delta_k(k*self.h, redshift)/k**3

    def sigma_r(self, scale, redshift=None):
        """RMS power on scale in Mpc/h"""
        sigma2, sigma2_err = integrate.quad(
            self._sigma_integrand, numpy.log10(self.k_min),
            numpy.log10(self.k_max), args=(scale,))
        sigma2 *= numpy.log(10.0)/(2.0*numpy.pi**2)

        if not redshift is None:
            sigma2 *= self.growth_factor(redshift)**2

        return numpy.sqrt(sigma2)

    def _sigma_integrand(self, logk, scale):
        k = 10**logk
        kR = scale*k

        W = 9.0*(numpy.sin(kR) - kR*numpy.cos(kR))**2/(kR**6)

        return self.linear_power(k, 0.0)*W*k**3

    def sigma_m(self, mass, redshift=None):
        """RMS power on scale subtended by total mean mass in solar masses."""
        z = redshift
        if z is None: z = 0.0
        scale = (3.0*mass/(4.0*numpy.pi*self.rho_bar(z)))**(1.0/3.0)

        return self.sigma_r(scale, redshift)

    def nu_r(self, scale, redshift=None):
        """Ratio of (delta_c(z)/sigma(R, z))^2"""
        return (self.delta_c(redshift)/self.sigma_r(scale, redshift))**2.0

    def nu_m(self, mass, redshift=None):
        """Ratio of (delta_c(z)/sigma(M, z))^2"""
        return (self.delta_c(redshift)/self.sigma_m(mass, redshift))**2.0

    def write(self, output_file_name, output_power_file_name=None):
        f = open(output_file_name, "w")
        for z, chi, growth in zip(
            self._z_array, self._chi_array, self._growth_array):
            if z < self.z_max:
                f.write("%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n" % (
                    z, chi, growth, self.omega_m(z), self.omega_l(z),
                    self.delta_v(z), self.delta_c(z)))
        f.close()

        if not output_power_file_name is None:
            k_min = 0.001
            k_max = 100.0

            dln_k = (numpy.log(k_max) - numpy.log(k_min))/200
            ln_k_max = numpy.log(k_max) + dln_k
            ln_k_min = numpy.log(k_min) - dln_k

            ln_k_array = numpy.arange(ln_k_min, ln_k_max + dln_k, dln_k)
            f = open(output_power_file_name, "w")
            for ln_k in ln_k_array:
                k = numpy.exp(ln_k)
                f.write("%1.10f %1.10f\n" % (k, self.linear_power(k)))
            f.close()

