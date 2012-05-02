from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import integrate
from scipy import special
import copy
import param  # cosmological parameter object from Cosmopy
import defaults
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

__author__ = ("Ryan Scranton <ryan.scranton@gmail.com>, "+
              "Chris Morrison <morrison.chrisb@gmail.com>")

class Halo(object):
    """Basic halo model object.

    Given an HOD object and cosmological and halo parameters, a halo object
    should be able to generate a non-linear power spectrum for dark matter,
    galaxies or their cross-spectrum.  It should also be able to return
    the halo profile as a function of mass and radius and its Fourier transform
    as a function of mass and wavenumber.
    """
    def __init__(self, input_hod=None, redshift=None, cosmo_dict=None,
                 halo_dict=None, use_camb=False, **kws):
        # Hard coded, but we shouldn't expect halos outside of this range.
        self._k_min = 0.001
        self._k_max = 100.0
        ln_mass_min = numpy.log(1.0e9)
        ln_mass_max = numpy.log(5.0e16)
        
        dln_k = (numpy.log(self._k_max) - numpy.log(self._k_min))/200.0
        self.ln_k_max = numpy.log(self._k_max) + dln_k
        self.ln_k_min = numpy.log(self._k_min) - dln_k

        self._ln_k_array = numpy.arange(
            self.ln_k_min, self.ln_k_max + dln_k, dln_k)

        if cosmo_dict is None:
            cosmo_dict = defaults.default_cosmo_dict
        self.cosmo_dict = cosmo_dict

        if halo_dict is None:
            halo_dict = defaults.default_halo_dict
        self.halo_dict = halo_dict

        if redshift is None: redshift = 0.0
        self.redshift = redshift

        self.c0 = halo_dict["c0"]/(1.0 + self.redshift)
        self.beta = halo_dict["beta"]

        # If we hard-code to an NFW profile, then we can use an analytic
        # form for the halo profile Fourier transform.
        # self.alpha = -1.0*halo_dict.dpalpha
        self.alpha = -1.0

        self.cosmo = cosmology.SingleEpoch(self.redshift, cosmo_dict)
        self.delta_v = self.cosmo.delta_v()
        self.rho_bar = self.cosmo.rho_bar()
        self.h = self.cosmo.h

        self.mass = mass_function.MassFunction(
            self.redshift, cosmo_dict, halo_dict)

        if input_hod is None:
            input_hod = hod.HODZheng()
        self.local_hod = input_hod
    
        self.use_camb = use_camb
        if self.use_camb:
            self.camb = camb.CambWrapper(cosmo_dict)
            self.camb.set_redshift(redshift)
            self.camb.run()
            self.camb.normalize(self.cosmo.sigma_8*self.cosmo._growth)

        self._calculate_n_bar()
        self._initialize_halo_splines()

        self._initialized_h_m = False
        self._initialized_h_g = False

        self._initialized_pp_mm = False
        self._initialized_pp_gm = False
        self._initialized_pp_gg = False

    def set_cosmology(self, cosmo_dict=None, redshift=None):
        if cosmo_dict==None:
            cosmo_dict = self.cosmo_dict
        if redshift==None:
            redshift = self.redshift
        self.cosmo = cosmology.SingleEpoch(redshift, cosmo_dict)
        self.delta_v = self.cosmo.delta_v()
        self.rho_bar = self.cosmo.rho_bar()
        self.h = self.cosmo.h

        self.c0 = self.halo_dict["c0"]/(1.0 + redshift)

        self.mass = mass_function.MassFunction(
            self.redshift, cosmo_dict, self.halo_dict)

        self._calculate_n_bar()
        self._initialize_halo_splines()

        self._initialized_h_m = False
        self._initialized_h_g = False

        self._initialized_pp_mm = False
        self._initialized_pp_gm = False
        self._initialized_pp_gg = False

        if self.use_camb:
            self.camb = camb.CambWrapper(cosmo_dict)
            self.camb.set_redshift(redshift)
            self.camb.run()
            self.camb.normalize(self.cosmo.sigma_8*self.cosmo._growth)

    def set_hod(self, input_hod):
        self.local_hod = input_hod

        self._calculate_n_bar()

        self._initialized_h_m = False
        self._initialized_h_g = False

        self._initialized_pp_mm = False
        self._initialized_pp_gm = False
        self._initialized_pp_gg = False

    def set_halo(self, halo_dict=None):
        self.c0 = halo_dict["c0"]/(1.0 + self.redshift)
        self.beta = halo_dict["beta"]
        self.alpha = -1.0

        self.mass = mass_function.MassFunction(
            self.redshift, self.cosmo_dict, halo_dict)

        self.local_hod.set_halo(halo_dict)
        self.set_hod(self.local_hod)
        
    def set_redshift(self, redshift):
        self.set_cosmology(self.cosmo_dict, redshift)

    def linear_power(self, k):
        if self.use_camb:
            return self.linear_power_camb(k)
        else:
            return self.linear_power_eh(k)

    def linear_power_camb(self, k):
        """Linear power spectrum in comoving (Mpc/h)^3 from CAMB."""
        return self.camb.linear_power(k)

    def linear_power_eh(self, k):
        """Linear power spectrum in comoving (Mpc/h)^3 from E+Hu99."""
        return self.cosmo.linear_power(k)

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
        if not self._initialized_pp_gg:
            self._initialize_pp_gg()

        return self.linear_power(k)*self._h_g(k)*self._h_g(k) + self._pp_gg(k)

    def virial_radius(self, mass):
        """Halo virial radius in Mpc as a function of mass in M_sun."""
        return numpy.exp(self._ln_r_v_spline(numpy.log(mass)))

    def concentration(self, mass):
        """Halo concentration as a function of mass in M_sun."""
        return numpy.exp(self._ln_concen_spline(numpy.log(mass)))

    def halo_normalization(self, mass):
        """Halo concentration in M_sun/Mpc^3 as a function of mass in M_sun."""
        return numpy.exp(self._ln_halo_norm_spline(numpy.log(mass)))

    def y(self, ln_k, mass):
        """Fourier transform of the halo profile.

        Using an analytic expression for the Fourier transform from
        White (2001).  This is only valid for an NFW profile.
        """

        k = numpy.exp(ln_k)/self.h
        con = self.concentration(mass)
        con_plus = 1.0 + con
        z = k*self.virial_radius(mass)/con
        si_z, ci_z = special.sici(z)
        si_cz, ci_cz = special.sici(con_plus*z)
        rho_km = (numpy.cos(z)*(ci_cz - ci_z) +
                  numpy.sin(z)*(si_cz - si_z) -
                  numpy.sin(con*z)/(con_plus*z))
        mass_k = numpy.log(con_plus) - con/con_plus

        return rho_km/mass_k

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
                numpy.log10(mass), self.y(ln_k, mass), self.concentration(mass),
                self.halo_normalization(mass), self.virial_radius(mass)))
        f.close()

    def write_power_components(self, output_file_name):
        f = open(output_file_name, "w")
        for ln_k in self._ln_k_array:
            k = numpy.exp(ln_k)
            f.write("%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n" % (
                k, self._h_m(k), self._pp_mm(k),
                self._h_g(k), self._pp_gm(k), self._pp_gg(k)))
        f.close()

    def _h_m(self, k):
        if k >= self._k_min and k <= self._k_max:
            return self._h_m_spline(numpy.log(k))
        else:
            return 0.0

    def _pp_mm(self, k):
        if k >= self._k_min and k <= self._k_max:
            return self._pp_mm_spline(numpy.log(k))
        else:
            return 0.0

    def _pp_gm(self, k):
        if k >= self._k_min and k <= self._k_max:
            return self._pp_gm_spline(numpy.log(k))
        else:
            return 0.0

    def _h_g(self, k):
        if k >= self._k_min and k <= self._k_max:
            return self._h_g_spline(numpy.log(k))
        else:
            return 0.0

    def _pp_gg(self, k):
        if k >= self._k_min and k <= self._k_max:
            return self._pp_gg_spline(numpy.log(k))
        else:
            return 0.0

    def _calculate_n_bar(self):
        self.n_bar_over_rho_bar, nbar_err = integrate.quad(
            self._nbar_integrand, numpy.log(self.mass.nu_min),
            numpy.log(self.mass.nu_max),limit=200)

        self.n_bar = self.n_bar_over_rho_bar*self.rho_bar

    def _nbar_integrand(self, ln_nu):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)

        return nu*self.local_hod.first_moment(mass)*self.mass.f_nu(nu)/mass

    # def _calculate_n_bar(self):
    #     self.n_bar, nbar_err = integrate.quad(
    #         self._nbar_integrand, self.mass.ln_mass_min,
    #         self.mass.ln_mass_max)

    #     self.n_bar_over_rho_bar = self.n_bar/self.rho_bar

    # def _nbar_integrand(self, ln_mass):
    #     mass = self.mass.mass(ln_mass)

    #     return mass*self.local_hod.first_moment(mass)*self.mass.f_m(mass)

    def _calculate_bias(self):
        self.bias, bias_err = integrate.quad(
            self._bias_integrand, self.mass.ln_mass_min,self.mass.ln_mass_max,
            limit=200)

        self.bias = self.bias/self.n_bar

    def _bias_integrand(self, ln_mass):
        mass = numpy.exp(ln_mass)

        return (mass*self.local_hod.first_moment(mass)*self.mass.f_m(mass)*
                self.mass.bias_m(mass))

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

    def a_c(self, mass):
        """Formation epoch definition from Wechsler et al. 2002
        """
        a_c = 0.1*numpy.log10(mass)-0.9
        if a_c > 0.01:
            return 0.1*numpy.log10(mass)-0.9
        elif a_c <=0.01:
            return 0.01

    # def _concentration(self, mass):
    #     """Halo concentration as a function of halo mass.

    #     Functional form from Wechsler et al. 2002
    #     """
    #     return 4.1/(self.a_c(mass)*(1+self.redshift))

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
                numpy.log(self.mass.nu_max),limit=200,
                args=(self._ln_k_array[idx],))
            h_m_array[idx] = h_m

        self._h_m_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, h_m_array)
        self._initialized_h_m = True

    def _h_m_integrand(self, ln_nu, ln_k):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)

        return nu*self.mass.f_nu(nu)*self.mass.bias_nu(nu)*self.y(ln_k, mass)

    def _initialize_h_g(self):
        h_g_array = numpy.zeros_like(self._ln_k_array)

        for idx in xrange(self._ln_k_array.size):
            h_g, h_g_err = integrate.quad(
                self._h_g_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=200,
                args=(self._ln_k_array[idx],)) #Limit Was 150
            h_g_array[idx] = h_g/self.n_bar_over_rho_bar

        self._h_g_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, h_g_array)
        self._initialized_h_g = True

    def _h_g_integrand(self, ln_nu, ln_k):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)

        return (nu*self.mass.f_nu(nu)*self.mass.bias_nu(nu)*
                self.y(ln_k, mass)*self.local_hod.first_moment(mass)/mass)

    def _initialize_pp_mm(self):
        pp_mm_array = numpy.zeros_like(self._ln_k_array)

        for idx in xrange(self._ln_k_array.size):
            pp_mm, pp_mm_err = integrate.quad(
                self._pp_mm_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=200,
                args=(self._ln_k_array[idx],)) # Limit was 150
            pp_mm_array[idx] = pp_mm/self.rho_bar

        self._pp_mm_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, pp_mm_array)
        self._initialized_pp_mm = True

    def _pp_mm_integrand(self, ln_nu, ln_k):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)
        y = self.y(ln_k, mass)

        return nu*self.mass.f_nu(nu)*mass*y*y

    def _initialize_pp_gg(self):
        pp_gg_array = numpy.zeros_like(self._ln_k_array)

        for idx in xrange(self._ln_k_array.size):
            pp_gg, pp_gg_err = integrate.quad(
                self._pp_gg_integrand, numpy.log(self.mass.nu_min),
                numpy.log(self.mass.nu_max),limit=200,
                args=(self._ln_k_array[idx],))
            pp_gg_array[idx] = pp_gg*self.rho_bar/(self.n_bar*self.n_bar)

        self._pp_gg_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, pp_gg_array)
        self._initialized_pp_gg = True

    def _pp_gg_integrand(self, ln_nu, ln_k):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)
        y = self.y(ln_k, mass)
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
                numpy.log(self.mass.nu_max),limit=200,
                args=(self._ln_k_array[idx],)) # Limit Was 150
            pp_gm_array[idx] = pp_gm/self.n_bar

        self._pp_gm_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, pp_gm_array)
        self._initialized_pp_gm = True

    def _pp_gm_integrand(self, ln_nu, ln_k):
        nu = numpy.exp(ln_nu)
        mass = self.mass.mass(nu)
        y = self.y(ln_k, mass)
        n_exp = self.local_hod.first_moment(self.mass.mass(nu))

        if n_exp < 1:
            return nu*self.mass.f_nu(nu)*n_exp*y
        else:
            return nu*self.mass.f_nu(nu)*n_exp*y*y


class HaloExclusion(Halo):

    def __init__(self, input_hod=None, redshift=None, cosmo_dict=None,
                 halo_dict=None, **kws):
        Halo.__init__(self, input_hod, redshift, cosmo_dict, halo_dict, **kws)
        ln_r_v_array = numpy.zeros_like(self.mass._nu_array)

        for idx in xrange(self.mass._nu_array.size):
            mass = numpy.exp(self.mass._ln_mass_array[idx])
            ln_r_v_array[idx] = numpy.log(self._virial_radius(mass))

        self._ln_nu_v_r_spline = InterpolatedUnivariateSpline(
            ln_r_v_array, self.mass._nu_array)

        self.v_r_max = numpy.exp(numpy.max(ln_r_v_array))
        self.v_r_min = numpy.exp(numpy.min(ln_r_v_array))
        self._initialized_h_m_ext = False

    def power_gm(self, k):
        """Galaxy-matter cross-spectrum in comoving (Mpc/h)^3"""
        if not self._initialized_h_m_ext:
            self._initialize_h_m_ext()
        if not self._initialized_h_g:
            self._initialize_h_g()
        if not self._initialized_pp_gm:
            self._initialize_pp_gm()

        return (self.power_mm(k)*self._h_g(k)*self._h_m_ext(k) + 
                self._pp_gm(k))

    def power_mg(self, k):
        """Galaxy-matter cross-spectrum in comoving (Mpc/h)^3"""
        return self.power_gm(k)

    def power_gg(self, k):
        """Galaxy power spectrum in comoving (Mpc/h)^3"""
        if not self._initialized_h_g:
            self._initialize_h_g()
        if not self._initialized_pp_gg:
            self._initialize_pp_gg()

        return (self.power_mm(k)*self._h_g(k)*self._h_g(k) + 
                self._pp_gg(k))

    def _h_m_ext(self, k):
        if k >= self._k_min and k <= self._k_max:
            return self._h_m_ext_spline(numpy.log(k))
        else:
            return 0.0

    def _initialize_h_m_ext(self):
        h_m_ext_array = numpy.zeros_like(self._ln_k_array)
        
        for idx in xrange(self._ln_k_array.size):
            r_lim = 0.5/numpy.exp(self._ln_k_array[idx])
            if r_lim >= self.v_r_max:
                r_lim = self.v_r_max
            elif r_lim <= self.v_r_min:
                r_lim = self.v_r_min
            nu_max = self._ln_nu_v_r_spline(numpy.log(r_lim))
            if nu_max >= self.mass.nu_max:
                nu_max = self.mass.nu_max
            elif nu_max <= self.mass.nu_min:
                nu_max = self.mass.nu_min
            
            h_m, h_m_err = integrate.quad(
                self._h_m_integrand, numpy.log(self.mass.nu_min),
                numpy.log(nu_max),limit=200,
                args=(self._ln_k_array[idx],))
            h_m_ext_array[idx] = h_m

        self._h_m_ext_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, h_m_ext_array)
        self._initialized_h_m_ext = True

    def _initialize_h_g(self):
        h_g_array = numpy.zeros_like(self._ln_k_array)

        for idx in xrange(self._ln_k_array.size):
            r_lim = 0.5/numpy.exp(self._ln_k_array[idx])
            if r_lim >= self.v_r_max:
                r_lim = self.v_r_max
            elif r_lim <= self.v_r_min:
                r_lim = self.v_r_min
            nu_max = self._ln_nu_v_r_spline(numpy.log(r_lim))
            if nu_max >= self.mass.nu_max:
                nu_max = self.mass.nu_max
            elif nu_max <= self.mass.nu_min:
                nu_max = self.mass.nu_min

            h_g, h_g_err = integrate.quad(
                self._h_g_integrand, numpy.log(self.mass.nu_min),
                numpy.log(nu_max),limit=200,
                args=(self._ln_k_array[idx],))
            h_g_array[idx] = h_g/self.n_bar_over_rho_bar

        self._h_g_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, h_g_array)
        self._initialized_h_g = True
