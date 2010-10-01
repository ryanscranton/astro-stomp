from scipy.interpolate import SmoothBivariateSpline
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
        self._k_min = 0.001
        self._k_max = 100.0

        dln_k = (numpy.log(self._k_max) - numpy.log(self._k_min))/200.0
        self.ln_k_max = numpy.log(self._k_max) + dln_k
        self.ln_k_min = numpy.log(self._k_min) - dln_k

        self._ln_k_array = numpy.arange(
            self.ln_k_min, self.ln_k_max + dln_k, dln_k)

        if camb_param is None:
            camb_param = param.CambParams(**kws)
        self.camb_param = camb_param

        if halo_param is None:
            halo_param = param.HaloModelParams(**kws)
        self.halo_param = halo_param

        if redshift is None: redshift = 0.0
        self.redshift = redshift

        self.c0 = halo_param.cbarcoef/(1.0 + self.redshift)
        self.beta = halo_param.cbarslope
        self.alpha = -1.0*halo_param.dpalpha

        cosmo = cosmology.SingleEpoch(self.redshift, camb_param)
        self.delta_v = cosmo.delta_v()
        self.rho_bar = cosmo.rho_bar()
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
        self._initialize_halo_splines()
        self._initialize_y_splines()

        self._initialized_h_m = False
        self._initialized_h_g = False

        self._initialized_pp_mm = False
        self._initialized_pp_gm = False
        self._initialized_pp_gg = False

    def set_cosmology(self, camb_param=None, redshift=None):
        if camb_param==None:
            camb_param = self.camb_param
        if redshift==None:
            redshift = self.redshift
        cosmo = cosmology.SingleEpoch(self.redshift, camb_param)
        self.delta_v = cosmo.delta_v()
        self.rho_bar = cosmo.rho_bar()
        self.h = cosmo.h

        self.c0 = self.halo_param.cbarcoef/(1.0 + self.redshift)

        self.mass = mass_function.MassFunction(
            self.redshift, camb_param, self.halo_param)

        self.camb = camb.CambWrapper(camb_param)
        self.camb.set_redshift(self.redshift)
        self.camb.run()

        self._calculate_n_bar()
        self._initialize_halo_splines()
        self._initialize_y_splines()

        self._initialized_h_m = False
        self._initialized_h_g = False

        self._initialized_pp_mm = False
        self._initialized_pp_gm = False
        self._initialized_pp_gg = False

    def set_hod(self, input_hod):
        self.local_hod = input_hod

        self._calculate_n_bar()
        self._initialize_halo_splines()
        self._initialize_y_splines()

        self._initialized_h_m = False
        self._initialized_h_g = False

        self._initialized_pp_mm = False
        self._initialized_pp_gm = False
        self._initialized_pp_gg = False

    def set_halo(self, halo_param=None):
        self.c0 = halo_param.cbarcoef/(1.0 + self.redshift)
        self.beta = halo_param.cbarslope
        self.alpha = -1.0*halo_param.dpalpha

        self.mass = mass_function.MassFunction(
            self.redshift, self.camb_param, halo_param)

        self.local_hod.set_halo(halo_param)
        self.set_hod(self.local_hod)
        
    def set_redshift(self, redshift):
        self.set_cosmology(self.camb_param, redshift)

    def linear_power(self, k):
        """Linear power spectrum in comoving (Mpc/h)^3 from CAMB."""
        return self.camb.linear_power(k)

    def power_mm(self, k):
        """Non-linear power spectrum in comoving (Mpc/h)^3"""
        if not self._initialized_h_m:
            self._initialize_h_m()
        if not self._initialized_pp_mm:
            self._initialize_pp_mm()

        return self.linear_power(k)*self._h_m(k)*self._h_m(k) + self._pp_mm(k)

    def power_gm(self, k):
        """Galaxy-matter cross-spectrum in comoving (Mpc/h)^3"""
        if not self._initialized_h_m:
            self._initialize_h_m()
        if not self._initialized_h_g:
            self._initialize_h_g()
        if not self._initialized_pp_gm:
            self._initialize_pp_gm()

        return self.linear_power(k)*self._h_g(k)*self._h_m(k) + self._pp_gm(k)

    def power_mg(self, k):
        """Galaxy-matter cross-spectrum in comoving (Mpc/h)^3"""
        return self.power_gm(k)

    def power_gg(self, k):
        """Galaxy power spectrum in comoving (Mpc/h)^3"""
        if not self._initialized_h_g:
            self._initialize_h_g()
        if not self._initialized_pp_gm:
            self._initialize_pp_gg()

        return self.linear_power(k)*self._h_g(k)*self._h_g(k) + self._pp_gg(k)

    def virial_radius(self, mass):
        """Halo virial radius in Mpc as a function of mass in M_sun."""
        return numpy.exp(self._ln_r_v_spline(numpy.log(mass))[0])

    def concentration(self, mass):
        """Halo concentration as a function of mass in M_sun."""
        return numpy.exp(self._ln_concen_spline(numpy.log(mass))[0])

    def halo_normalization(self, mass):
        """Halo concentration in M_sun/Mpc^3 as a function of mass in M_sun."""
        return numpy.exp(self._ln_halo_norm_spline(numpy.log(mass))[0])

    def write(self, output_file_name):
        f = open(output_file_name, "w")
        for ln_k in self._ln_k_array:
            k = numpy.exp(ln_k)
            f.write("%1.10f %1.10f %1.10f %1.10f %1.10f\n" % (
                k, self.linear_power(k), self.power_mm(k),
                self.power_gg(k), self.power_gm(k)))
        f.close()

    def write_halo(self, output_file_name, k=None):
        if k is None:
            k = 0.001
        ln_k = numpy.log(k)

        f = open(output_file_name, "w")
        for nu in self.mass._nu_array:
            mass = self.mass.mass(nu)
            f.write("%1.10f %1.10f %1.10f %1.10f %1.10f\n" % (
                numpy.log10(mass), self._y(ln_k, nu), self.concentration(mass),
                self.halo_normalization(mass), self.virial_radius(mass)))
        f.close()

    def _h_m(self, k):
        if k >= self._k_min and k <= self._k_max:
            return self._h_m_spline(numpy.log(k))[0]
        else:
            return 0.0

    def _pp_mm(self, k):
        if k >= self._k_min and k <= self._k_max:
            return self._pp_gm_spline(numpy.log(k))[0]
        else:
            return 0.0

    def _pp_gm(self, k):
        if k >= self._k_min and k <= self._k_max:
            return self._pp_gm_spline(numpy.log(k))[0]
        else:
            return 0.0

    def _h_g(self, k):
        if k >= self._k_min and k <= self._k_max:
            return self._h_g_spline(numpy.log(k))[0]
        else:
            return 0.0

    def _pp_gg(self, k):
        if k >= self._k_min and k <= self._k_max:
            return self._pp_gg_spline(numpy.log(k))[0]
        else:
            return 0.0

    def _calculate_n_bar(self):
        self.n_bar_over_rho_bar, nbar_err = integrate.quad(
            self._nbar_integrand, numpy.log(self.mass.nu_min),
            numpy.log(self.mass.nu_max))

        self.n_bar = self.n_bar_over_rho_bar*self.rho_bar

    def _nbar_integrand(self, ln_nu):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)

        return nu*self.local_hod.first_moment(mass)*self.mass.f_nu(nu)/mass

    def _initialize_halo_splines(self):
        ln_r_v_array = numpy.zeros_like(self.mass._nu_array)
        ln_concen_array = numpy.zeros_like(self.mass._nu_array)
        ln_halo_norm_array = numpy.zeros_like(self.mass._nu_array)

        for idx in xrange(self.mass._nu_array.size):
            mass = numpy.exp(self.mass._ln_mass_array[idx])
            ln_r_v_array[idx] = numpy.log(self._virial_radius(mass))
            ln_concen_array[idx] = numpy.log(self._concentration(mass))
            ln_halo_norm_array[idx] = numpy.log(self._halo_normalization(mass))

        self._ln_r_v_spline = InterpolatedUnivariateSpline(
            self.mass._ln_mass_array, ln_r_v_array)
        self._ln_concen_spline = InterpolatedUnivariateSpline(
            self.mass._ln_mass_array, ln_concen_array)
        self._ln_halo_norm_spline = InterpolatedUnivariateSpline(
            self.mass._ln_mass_array, ln_halo_norm_array)

    def _concentration(self, mass):
        """Halo concentration as a function of halo mass.

        Functional form from Bullock et al.
        """
        return self.c0*(mass/self.mass.m_star)**self.beta

    def _halo_normalization(self, mass):
        """Halo normalization as a function of mass.

        The halo density profile is normalized such that the integral of the
        density profile out to the virial radius equals the halo mass.  This
        ends up being the ratio between the halo mass and the integral

        int(0, concentration, x**(2+alpha)/(1 + x)**(3+alpha))

        which is a hypergeometric function of alpha and the concentration.
        """
        con = self._concentration(mass)
        rho_s = (self.rho_bar*self.delta_v*con*con*con)/3.0
        rho_norm = (con**(3.0 + self.alpha)*
                    special.hyp2f1(3.0+self.alpha, 3.0+self.alpha,
                                   4.0+self.alpha, -1.0*con))/(3.0+self.alpha)

        return rho_s/rho_norm

    def _virial_radius(self, mass):
        """Halo virial radius as a function of mass."""
        r3 = 3.0*mass/(4.0*numpy.pi*self.delta_v*self.rho_bar)
        return r3**(1.0/3.0)

    def _initialize_h_m(self):
        h_m_array = numpy.zeros_like(self._ln_k_array)

        for idx in xrange(self._ln_k_array.size):
            h_m, h_m_err = integrate.quad(
                self._h_m_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=150,
                args=(self._ln_k_array[idx],))
            h_m_array[idx] = h_m

        self._h_m_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, h_m_array)
        self._initialized_h_m = True

    def _h_m_integrand(self, ln_nu, ln_k):
        nu = numpy.exp(ln_nu)

        return nu*self.mass.f_nu(nu)*self.mass.bias_nu(nu)*self._y(ln_k, nu)

    def _initialize_h_g(self):
        h_g_array = numpy.zeros_like(self._ln_k_array)

        for idx in xrange(self._ln_k_array.size):
            h_g, h_g_err = integrate.quad(
                self._h_g_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=150,
                args=(self._ln_k_array[idx],))
            h_g_array[idx] = h_g/self.n_bar_over_rho_bar

        self._h_g_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, h_g_array)
        self._initialized_h_g = True

    def _h_g_integrand(self, ln_nu, ln_k):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)

        return (nu*self.mass.f_nu(nu)*self.mass.bias_nu(nu)*
                self._y(ln_k, nu)*self.local_hod.first_moment(mass)/mass)

    def _initialize_pp_mm(self):
        pp_mm_array = numpy.zeros_like(self._ln_k_array)

        for idx in xrange(self._ln_k_array.size):
            pp_mm, pp_mm_err = integrate.quad(
                self._pp_mm_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=150,
                args=(self._ln_k_array[idx],))
            pp_mm_array[idx] = pp_mm/self.rho_bar

        self._pp_mm_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, pp_mm_array)
        self._initialized_pp_mm = True

    def _pp_mm_integrand(self, ln_nu, ln_k):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)
        y = self._y(ln_k, nu)

        return nu*self.mass.f_nu(nu)*mass*y*y

    def _initialize_pp_gg(self):
        pp_gg_array = numpy.zeros_like(self._ln_k_array)

        for idx in xrange(self._ln_k_array.size):
            pp_gg, pp_gg_err = integrate.quad(
                self._pp_gg_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=150,
                args=(self._ln_k_array[idx],))
            pp_gg_array[idx] = pp_gg*self.rho_bar/(self.n_bar*self.n_bar)

        self._pp_gg_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, pp_gg_array)
        self._initialized_pp_gg = True

    def _pp_gg_integrand(self, ln_nu, ln_k):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)
        y = self._y(ln_k, nu)
        n_pair = self.local_hod.second_moment(self.mass.mass(nu))

        if n_pair < 1:
            return nu*self.mass.f_nu(nu)*n_pair*y/mass
        else:
            return nu*self.mass.f_nu(nu)*n_pair*y*y/mass

    def _initialize_pp_gm(self):
        pp_gm_array = numpy.zeros_like(self._ln_k_array)

        for idx in xrange(self._ln_k_array.size):
            pp_gm, pp_gm_err = integrate.quad(
                self._pp_gm_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=150,
                args=(self._ln_k_array[idx],))
            pp_gm_array[idx] = pp_gm/self.n_bar

        self._pp_gm_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, pp_gm_array)
        self._initialized_pp_gm = True

    def _pp_gm_integrand(self, ln_nu, ln_k):
        nu = numpy.exp(ln_nu)
        y = self._y(ln_k, nu)
        n_exp = self.local_hod.first_moment(self.mass.mass(nu))

        if n_exp < 1:
            return nu*self.mass.f_nu(nu)*n_exp*y
        else:
            return nu*self.mass.f_nu(nu)*n_exp*y*y

    def _y(self, ln_k, nu):
        return self._y_spline(ln_k, nu)[0]

    def _initialize_y_splines(self):
        """Calculate y(k, mass) for a given wavenumber and mass.

        This is a time intensive part of the code, but we can speed things up
        by truncating the calculations when y -> 0 as k >> 1 and M >> M_star.
        """

        dln_k = (numpy.log(self._k_max) - numpy.log(self._k_min))/50
        ln_k_array = numpy.arange(
            self.ln_k_min, self.ln_k_max + dln_k, dln_k)

        dln_mass = (self.mass.ln_mass_max - self.mass.ln_mass_min)/50
        ln_mass_array = numpy.arange(
            self.mass.ln_mass_min, self.mass.ln_mass_max + dln_mass, dln_mass)

        self.y_array = numpy.zeros(ln_k_array.size*ln_mass_array.size)
        self.y_k_array = numpy.zeros(ln_k_array.size*ln_mass_array.size)
        self.y_nu_array = numpy.zeros(ln_k_array.size*ln_mass_array.size)

        tol = 0.001

        kdx = 0
        for idx in xrange(ln_k_array.size):
            k = numpy.exp(ln_k_array[idx])

            calculate_y = True
            for jdx in xrange(ln_mass_array.size):
                nu = self.mass.nu(numpy.exp(ln_mass_array[jdx]))
                if calculate_y:
                    mass = numpy.exp(ln_mass_array[jdx])
                    r_v = self.virial_radius(mass)
                    y, y_err = integrate.quad(
                        self._y_integrand,
                        numpy.log(1.0e-6*r_v), numpy.log(r_v),
                        args=(k, mass,), limit=150)
                    self.y_array[kdx] = y/mass
                else:
                    self.y_array[kdx] = 0.0

                if self.y_array[kdx] < tol:
                    calculate_y = False
                self.y_k_array[kdx] = ln_k_array[idx]
                self.y_nu_array[kdx] = nu
                kdx += 1

        self._y_spline = SmoothBivariateSpline(self.y_k_array, self.y_nu_array, self.y_array)

    def _y_integrand(self, ln_radius, k, mass):
        r = numpy.exp(ln_radius)
        k /= self.h
        r_norm = (r*self.concentration(mass)/self.virial_radius(mass))
        rho = r_norm**self.alpha/(1.0 + r_norm)**(3.0 + self.alpha)
        rho *= self.halo_normalization(mass)

        return 4.0*numpy.pi*r*r*r*rho*numpy.sin(k*r)/(k*r)
