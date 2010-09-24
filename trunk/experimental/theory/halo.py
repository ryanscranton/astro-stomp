from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import integrate
from scipy import special
import copy
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

__author__ = ("Ryan Scranton <ryan.scranton@gmail.com>, "
              "Chris Morrison <morrison.chrisb@gmail.com>")

class Halo(object):
    """Basic halo model object.

    Given an HOD object and cosmological and halo parameters, a halo object
    should be able to generate a non-linear power spectrum for dark matter,
    galaxies or their cross-spectrum.  It should also be able to return
    the halo profile as a function of mass and radius and its Fourier transform
    as a function of mass and wavenumber.
    """
    def __init__(self, input_hod=None, redshift=None, camb_param=None,
                 halo_param=None, **kws):
        # Hard coded, but we shouldn't expect halos outside of this range.
        k_min = 0.001
        k_max = 100.0

        dln_k = (numpy.log(k_max) - numpy.log(k_min))/200
        self.ln_k_max = numpy.log(k_max) + dln_k
        self.ln_k_min = numpy.log(k_min) - dln_k

        self._ln_k_array = numpy.arange(
            self.ln_k_min, self.ln_k_max + dln_k, dln_k)

        if camb_param is None:
            camb_param = param.CambParams(**kws)

        if halo_param is None:
            halo_param = param.HaloModelParams(**kws)

        if redshift is None: redshift = 0.0
        self.redshift = redshift
        z_min = 0.0
        z_max = self.redshift + 1.0

        self.c0 = halo_param.cbarcoef/(1.0 + self.redshift)
        self.beta = halo_param.cbarslope
        self.alpha = -1.0*halo_param.dpalpha

        cosmo = cosmology.Cosmology(z_min, z_max, camb_param)
        self.delta_v = cosmo.delta_v(self.redshift)
        self.rho_bar = cosmo.rho_bar(self.redshift)
        self.h = cosmo.h

        self.mass = mass_function.MassFunction(
            self.redshift, camb_param, halo_param)

        if input_hod is None:
            self.local_hod = hod.HODKravtsov(halo_param)
        else:
            self.local_hod = input_hod

        self.camb = camb.CambWrapper(camb_param)
        self.camb.set_redshift(redshift)
        self.camb.run()

        self._calculate_n_bar()
        self._initialize_splines()

    def linear_power(self, k):
        return self.camb.linear_power(k)

    def power_mm(self, k):
        return self._Phh_mm(k) + self._pp_mm(k)

    def power_gm(self, k):
        return self._Phh_gm(k) + self._pp_gm(k)

    def power_gg(self, k):
        return self._Phh_gg(k) + self._pp_gg(k)

    def _Phh_mm(self, k):
        return self.linear_power(k)*(self._h_m_spline(numpy.log(k)))**2

    def _pp_mm(self, k):
        return self._pp_gm_spline(numpy.log(k))[0]

    def _Phh_gm(self, k):
        return (self.linear_power(k)*self._h_m_spline(numpy.log(k))[0]*
                self._h_g_spline(numpy.log(k))[0])

    def _pp_gm(self, k):
        return self._pp_gm_spline(numpy.log(k))[0]

    def _Phh_gg(self, k):
        return self.linear_power(k)*(self._h_g_spline(numpy.log(k))[0])**2

    def _pp_gg(self, k):
        return self._pp_gg_spline(numpy.log(k))[0]

    def _calculate_n_bar(self):
        self.n_bar_over_rho_bar, nbar_err = integrate.quad(
            self._nbar_integrand, numpy.log(self.mass.nu_min),
            numpy.log(self.mass.nu_max))

        self.n_bar = self.n_bar_over_rho_bar*self.rho_bar

    def _nbar_integrand(self, ln_nu):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)

        return nu*self.local_hod.first_moment(mass)*self.mass.f_nu(nu)/mass

    def _initialize_splines(self):
        h_m_array = numpy.zeros_like(self._ln_k_array)
        h_g_array = numpy.zeros_like(self._ln_k_array)

        pp_mm_array = numpy.zeros_like(self._ln_k_array)
        pp_gg_array = numpy.zeros_like(self._ln_k_array)
        pp_gm_array = numpy.zeros_like(self._ln_k_array)

        for idx in xrange(self._ln_k_array.size):
            self._initialize_y(self._ln_k_array[idx])

            h_m, h_m_err = integrate.quad(
                self._h_m_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=150)
            h_m_array[idx] = h_m

            h_g, h_g_err = integrate.quad(
                self._h_g_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=150)
            h_g_array[idx] = h_g/self.n_bar_over_rho_bar

            pp_mm, pp_mm_err = integrate.quad(
                self._pp_mm_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=150)
            pp_mm_array[idx] = pp_mm/self.rho_bar

            pp_gg, pp_gg_err = integrate.quad(
                self._pp_gg_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=150)
            pp_gg_array[idx] = pp_gg*self.rho_bar/(self.n_bar*self.n_bar)

            pp_gm, pp_gm_err = integrate.quad(
                self._pp_gm_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=150)
            pp_gm_array[idx] = pp_gm/self.n_bar

        self._h_m_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, h_m_array)
        self._h_g_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, h_g_array)

        self._pp_mm_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, pp_mm_array)
        self._pp_gg_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, pp_gg_array)
        self._pp_gm_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, pp_gm_array)

    def _h_m_integrand(self, ln_nu):
        nu = numpy.exp(ln_nu)

        return nu*self.mass.f_nu(nu)*self.mass.bias_nu(nu)*self.y(nu)

    def _h_g_integrand(self, ln_nu):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)

        return (nu*self.mass.f_nu(nu)*self.mass.bias_nu(nu)*
                self.y(nu)*self.local_hod.first_moment(mass)/mass)

    def _pp_mm_integrand(self, ln_nu):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)

        return nu*self.mass.f_nu(nu)*mass*self.y(nu)*self.y(nu)

    def _pp_gg_integrand(self, ln_nu):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)
        y = self.y(nu)
        n_pair = self.local_hod.second_moment(self.mass.mass(nu))

        if n_pair < 1:
            return nu*self.mass.f_nu(nu)*n_pair*numpy.abs(y)/mass
        else:
            return nu*self.mass.f_nu(nu)*n_pair*y*y/mass

    def _pp_gm_integrand(self, ln_nu):
        nu = numpy.exp(ln_nu)
        y = self.y(nu)
        n_exp = self.local_hod.first_moment(self.mass.mass(nu))

        if n_exp < 1:
            return nu*self.mass.f_nu(nu)*n_exp*numpy.abs(y)
        else:
            return nu*self.mass.f_nu(nu)*n_exp*y*y

    def y(self, nu):
        return self._y_spline(nu)

    def _initialize_y(self, ln_k):
        y_array = numpy.zeros_like(self.mass._nu_array)

        for idx in xrange(y_array.size):
            mass = numpy.exp(self.mass._ln_mass_array[idx])
            r_v = self.virial_radius(mass)
            y, y_err = integrate.quad(self._y_integrand,
                                      numpy.log(1.0e-6*r_v), numpy.log(r_v),
                                      args=(numpy.exp(ln_k), mass,),
                                      limit=150)
            y_array[idx] = y/mass
        self._y_spline = InterpolatedUnivariateSpline(
            self.mass._nu_array, y_array)

    def _y_integrand(self, ln_radius, k, mass):
        radius = numpy.exp(ln_radius)
        k /= self.h
        r_norm = (radius*self.concentration(mass)/self.virial_radius(mass))
        rho = r_norm**self.alpha/(1.0 + r_norm)**(3.0+self.alpha)
        rho *= self.halo_normalization(mass)

        return 4.0*numpy.pi*radius*radius*rho*numpy.sin(k*radius)/(k*radius)

    def concentration(self, mass):
        """Halo concentration as a function of halo mass.

        Functional form from Bullock et al.
        """
        return self.c0*(mass/self.mass.m_star)**self.beta

    def halo_normalization(self, mass):
        """Halo normalization as a function of mass.

        The halo density profile is normalized such that the integral of the
        density profile out to the virial radius equals the halo mass.  This
        ends up being the ratio between the halo mass and the integral

        int(0, concentration, x**(2+alpha)/(1 + x)**(3+alpha))

        which is a hypergeometric function of alpha and the concentration.
        """
        con = self.concentration(mass)
        rho_s = (self.rho_bar*self.delta_v*con*con*con)/3.0
        rho_norm = (con**(3.0+self.alpha)*
                    special.hyp2f1(3.0+self.alpha, 3.0+self.alpha,
                                   4.0+self.alpha, -1.0*con))/(3.0+self.alpha)

        return rho_s/rho_norm

    def virial_radius(self, mass):
        """Halo virial radius as a function of mass."""
        r3 = 3.0*mass/(4.0*numpy.pi*self.delta_v*self.rho_bar)
        return r3**(1.0/3.0)

    def write(self, output_file_name):
        f = open(output_file_name, "w")
        for ln_k in self._ln_k_array:
            k = numpy.exp(ln_k)
            f.write("%1.10f %1.10f %1.10f %1.10f %1.10f\n" % (
                k, self.linear_power(k), self.power_mm(k),
                self.power_gg(k), self.power_gm(k)))
        f.close()

    def write_halo(self, output_file_name):
        f = open(output_file_name, "w")
        for nu in self.mass._nu_array:
            mass = self.mass.mass(nu)
            f.write("%1.10f %1.10f %1.10f %1.10f %1.10f\n" % (
                numpy.log10(mass), self.y(nu), self.concentration(mass),
                self.halo_normalization(mass), self.virial_radius(mass)))
        f.close()
