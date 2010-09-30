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


class MultiEpoch(object):
    """Class for calculating cosmological values across a range of redshifts.

    Given an initial set of redshift bounds and a set of parameters, MultiEpoch
    can return the comoving distance, angular diameter distance or luminosity
    distance for any redshift over the input interval.  Likewise, it can return
    a redshift, given an input comoving distance, as well as the growth factor
    D(z) and a number of other cosmological parameters.
    """
    def __init__(self, z_min, z_max, camb_param=None, **kws):
        dz = (z_max - z_min)/100.0
        self.z_max = z_max + dz
        self.z_min = z_min - dz
        if self.z_min < 0.0: self.z_min = 0.0

        dz = (self.z_max - self.z_min)/100.0
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

        self.flat = True
        self.open = False
        self.closed = False

        omega_total = self.omega_m0 + self.omega_l0
        if omega_total <= 1.0001 and omega_total >= 0.9999:
            self.flat = True
        else:
            self.flat = False
        if omega_total <= 0.9999:
            self.open = True
        else:
            self.open = False
        if omega_total > 1.0001:
            self.closed = True
        else:
            self.closed = False

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

    def comoving_distance(self, redshift):
        distance = 0.0

        if redshift < self.z_max and redshift > self.z_min:
            distance = self._chi_spline(redshift)[0]
        else:
            print ("Warning: requested redshift outside of bounds!  "
                   "Returning 0.")

        return distance

    def luminosity_distance(self, redshift):
        return (1.0 + redshift)*self.comoving_distance(redshift)

    def angular_diameter_distance(self):
        return self.comoving_distance(redshift)/(1.0 + redshift)

    def redshift(self, comoving_distance):
        return self._z_spline(comoving_distance)[0]

    def growth_factor(self, redshift):
        """Linear growth factor, normalized to unity at z = 0."""
        growth = 1.0

        if redshift < self.z_max and redshift > self.z_min:
            growth = self._growth_spline(redshift)[0]
        else:
            print ("Warning: requested redshift outside of bounds!  "
                   "Returning unity.")

        return growth

    def omega_m(self, redshift=None):
        if redshift is None: redshift = 0.0
        return self.omega_m0*(1.0 + redshift)**3/(
            self.omega_l0 + self.omega_r0*(1.0 + redshift)**2 +
            self.omega_m0*(1.0 + redshift)**3)

    def omega_l(self, redshift=None):
        if redshift is None: redshift = 0.0
        return self.omega_l0/(
            self.omega_l0 + self.omega_r0*(1.0 + redshift)**2 +
            self.omega_m0*(1.0 + redshift)**3)

    def delta_c(self, redshift=None):
        """Over-density threshold for linear, spherical collapse."""
        # Fitting function taken from NFW97
        delta_c = 0.15*(12.0*numpy.pi)**(2.0/3.0)
        if self.open:
            delta_c *= self.omega_m(redshift)**0.0185
        if self.flat and self.omega_m0 < 1.0001:
            delta_c *= self.omega_m(redshift)**0.0055

        if redshift is None:
            return delta_c
        else:
            return delta_c/self.growth_factor(redshift)

    def delta_v(self, redshift=None):
        """Over-density for a collapsed, virialized halo."""
        # Fitting function taken from NFW97
        delta_v = 178.0
        if self.open:
            delta_v /= self.omega_m(redshift)**0.7
        if self.flat and self.omega_m0 < 1.0001:
            delta_v /= self.omega_m(redshift)**0.55

        if redshift is None:
            return delta_v
        else:
            return delta_v/self.growth_factor(redshift)

    def rho_crit(self, redshift=None):
        """Critical density in solar masses per cubic Mpc."""
        if redshift is None: redshift = 0.0
        rho_crit = (self.omega_l0 + self.omega_r0*(1.0 + redshift)**2 +
                    self.omega_m0*(1.0 + redshift)**3)
        return rho_crit*1.0e-29*1.0e-33*2.937999e+73

    def rho_bar(self, redshift=None):
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
        return 2.0*numpy.pi*numpy.pi*self.delta_k(k*self.h, redshift)/(k*k*k)

    def sigma_r(self, scale, redshift=None):
        """RMS power on scale in Mpc/h"""
        sigma2, sigma2_err = integrate.quad(
            self._sigma_integrand, numpy.log(self.k_min),
            numpy.log(self.k_max), args=(scale,))
        sigma2 /= 2.0*numpy.pi*numpy.pi

        if not redshift is None:
            sigma2 *= self.growth_factor(redshift)**2

        return numpy.sqrt(sigma2)

    def _sigma_integrand(self, ln_k, scale):
        k = numpy.exp(ln_k)
        kR = scale*k

        W = 9.0*(numpy.sin(kR) - kR*numpy.cos(kR))**2/(kR*kR*kR*kR*kR*kR)

        return self.linear_power(k, 0.0)*W*k*k*k

    def sigma_m(self, mass, redshift=None):
        """RMS power on scale subtended by total mean mass in solar masses."""
        scale = (3.0*mass/(4.0*numpy.pi*self.rho_bar(redshift)))**(1.0/3.0)

        return self.sigma_r(scale, redshift)

    def nu_r(self, scale, redshift=None):
        """Ratio of (delta_c(z)/sigma(R, z))^2"""
        sqrt_nu = self.delta_c(redshift)/self.sigma_r(scale, redshift)
        return sqrt_nu*sqrt_nu

    def nu_m(self, mass, redshift=None):
        """Ratio of (delta_c(z)/sigma(M, z))^2"""
        sqrt_nu = self.delta_c(redshift)/self.sigma_m(mass, redshift)
        return sqrt_nu*sqrt_nu

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


class SingleEpoch(object):
    """Container class for calculating cosmological values at a given redshift.

    Given a redshift and a set of parameters, SingleEpoch can return the
    comoving distance, angular diameter distance or luminosity distance for.
    Likewise, it can return a redshift, given an input comoving distance, as
    well as the growth factor D(z) and other cosmological parameters.
    """
    def __init__(self, redshift, camb_param=None, **kws):
        self._redshift = redshift

        if camb_param is None:
            camb_param = param.CambParams(**kws)

        self.omega_m0 = camb_param.omega_cdm + camb_param.omega_baryon
        self.omega_b0 = camb_param.omega_baryon
        self.omega_l0 = camb_param.omega_lambda
        self.h = camb_param.hubble/100.0
        self.H0 = self.h/3000.0
        self.omega_r0 = 4.17e-5/(self.h**2)

        self.flat = True
        self.open = False
        self.closed = False

        omega_total = self.omega_m0 + self.omega_l0
        if omega_total <= 1.0001 and omega_total >= 0.9999:
            self.flat = True
        else:
            self.flat = False
        if omega_total <= 0.9999:
            self.open = True
        else:
            self.open = False
        if omega_total > 1.0001:
            self.closed = True
        else:
            self.closed = False

        self.k_min = 0.001
        self.k_max = 100.0
        self.n = camb_param.scalar_spectral_index[0]
        self.delta_H = (
            1.94e-5*self.omega_m0**(-0.785 - 0.05*numpy.log(self.omega_m0))*
            numpy.exp(-0.95*(self.n - 1) - 0.169*(self.n - 1)**2))

        self._initialize_defaults()

    def _initialize_defaults(self):
        self._chi, dist_err = integrate.quad(self.E, 0.0, self._redshift)

        a_min = 1.0e-5
        self.growth_norm, growth_norm_err = integrate.quad(
            self._growth_integrand, a_min, 1.0)
        self.growth_norm *= 2.5*self.omega_m0/(self.E(0.0)*self.H0)

        a = 1.0/(1.0 + self._redshift)
        growth, growth_err = integrate.quad(self._growth_integrand, a_min, a)
        growth *= 2.5*self.omega_m0/(self.E(self._redshift)*self.H0)
        self._growth = growth/self.growth_norm

    def E(self, redshift):
        return 1.0/(self.H0*numpy.sqrt(
            self.omega_l0 + self.omega_m0*(1 + redshift)**3 +
            (1 - self.omega_l0 - self.omega_m0)*(1 + redshift)**2))

    def _growth_integrand(self, redshift):
        return ((1.0 + redshift)*self.E(redshift)*self.H0)**3

    def comoving_distance(self):
        return self._chi

    def luminosity_distance(self):
        return (1.0 + self._redshift)*self._chi

    def angular_diameter_distance(self):
        return self._chi/(1.0 + self._redshift)

    def redshift(self):
        return self._redshift

    def growth_factor(self):
        """Linear growth factor, normalized to unity at z = 0."""
        return self._growth

    def omega_m(self):
        return self.omega_m0*(1.0 + self._redshift)**3/(
            self.omega_l0 + self.omega_r0*(1.0 + self._redshift)**2 +
            self.omega_m0*(1.0 + self._redshift)**3)

    def omega_l(self):
        return self.omega_l0/(
            self.omega_l0 + self.omega_r0*(1.0 + self._redshift)**2 +
            self.omega_m0*(1.0 + self._redshift)**3)

    def delta_c(self):
        """Over-density threshold for linear, spherical collapse."""
        # Fitting function taken from NFW97
        delta_c = 0.15*(12.0*numpy.pi)**(2.0/3.0)
        if self.open:
            delta_c *= self.omega_m()**0.0185
        if self.flat and self.omega_m0 < 1.0001:
            delta_c *= self.omega_m()**0.0055

        return delta_c/self._growth

    def delta_v(self):
        """Over-density for a collapsed, virialized halo."""
        # Fitting function taken from NFW97
        delta_v = 178.0
        if self.open:
            delta_v /= self.omega_m()**0.7
        if self.flat and self.omega_m0 < 1.0001:
            delta_v /= self.omega_m()**0.55

        return delta_v/self._growth

    def rho_crit(self):
        """Critical density in solar masses per cubic Mpc."""
        rho_crit = (self.omega_l0 + self.omega_r0*(1.0 + self._redshift)**2 +
                    self.omega_m0*(1.0 + self._redshift)**3)
        return rho_crit*1.0e-29*1.0e-33*2.937999e+73

    def rho_bar(self):
        """Matter density in solar masses per cubic Mpc."""
        return self.rho_crit()*self.omega_m()

    def _bbks(self, k):
        """BBKS transfer function."""
        q = k/self.omega_m0/self.h**2*numpy.exp(
            self.omega_b0 + numpy.sqrt(2*self.h)*self.omega_b0/self.omega_m0)
        return (numpy.log(1.0 + 2.34*q)/(2.34*q)*
                (1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**(-0.25))

    def delta_k(self, k):
        """k^2*P(k)/2*pi^2: dimensionless linear power spectrum."""
        delta_k = (
            self.delta_H**2*(3000.0*k/self.h)**(3 + self.n)*self._bbks(k)**2)
        return delta_k*self._growth*self._growth

    def linear_power(self, k):
        return 2.0*numpy.pi*numpy.pi*self.delta_k(k*self.h)/(k*k*k)

    def sigma_r(self, scale):
        """RMS power on scale in Mpc/h"""
        sigma2, sigma2_err = integrate.quad(
            self._sigma_integrand, numpy.log(self.k_min),
            numpy.log(self.k_max), args=(scale,), limit=100)
        sigma2 /= 2.0*numpy.pi*numpy.pi

        sigma2 *= self._growth*self._growth

        return numpy.sqrt(sigma2)

    def _sigma_integrand(self, ln_k, scale):
        k = numpy.exp(ln_k)
        kR = scale*k

        W = 9.0*(numpy.sin(kR) - kR*numpy.cos(kR))**2/(kR*kR*kR*kR*kR*kR)

        return self.linear_power(k)*W*k*k*k

    def sigma_m(self, mass):
        """RMS power on scale subtended by total mean mass in solar masses."""
        scale = (3.0*mass/(4.0*numpy.pi*self.rho_bar()))**(1.0/3.0)

        return self.sigma_r(scale)

    def nu_r(self, scale):
        """Ratio of (delta_c/sigma(R))^2"""
        sqrt_nu = self.delta_c()/self.sigma_r(scale)
        return sqrt_nu*sqrt_nu

    def nu_m(self, mass):
        """Ratio of (delta_c/sigma(M))^2"""
        sqrt_nu = self.delta_c()/self.sigma_m(mass)
        return sqrt_nu*sqrt_nu

    def write(self, output_power_file_name=None):
        print "z = %1.4f" % self._redshift
        print "Comoving distance = %1.4f" % self._chi
        print "Growth factor = %1.4f" % self._growth
        print "Omega_m(z) = %1.4f" % self.omega_m()
        print "Omega_l(z) = %1.4f" % self.omega_l()
        print "Delta_V(z) = %1.4f" % self.delta_v()
        print "delta_c(z) = %1.4f" % self.delta_c()
        print "sigma_8(z) = %1.4f" % self.sigma_r(8.0)

        if not output_power_file_name is None:
            dln_k = (numpy.log(self.k_max) - numpy.log(self.k_min))/200
            ln_k_max = numpy.log(self.k_max) + dln_k
            ln_k_min = numpy.log(self.k_min) - dln_k

            f = open(output_power_file_name, "w")
            for ln_k in numpy.arange(ln_k_min, ln_k_max + dln_k, dln_k):
                k = numpy.exp(ln_k)
                f.write("%1.10f %1.10f\n" % (k, self.linear_power(k)))
            f.close()

