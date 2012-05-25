import camb
import cosmology
import halo
import hod
import kernel
import defaults
import param
import numpy
from scipy import special
from scipy import integrate

"""Classes describing different correlation functions

Given a kernel function and halo model the class should produce a correlation as function of theta. The correlation should specifily the range of theta to be computed.

"""

speed_of_light = 3*10**5
degToRad = numpy.pi/180.0

__author__ = "Chris Morrison <morrison.chrisb@gmail.com>"

class Correlation(object):

    """
    Template class for correlation function, theta should be in radians. The integrand is over ln_k_theta. Choices of power spectra are 
    'power_gg' - The galaxy-galaxy power spectrum,
    'power_gm' - The galaxy-matter power spectrum,
    'power_mm' - The matter-matter power spectrum,
    'linear_power' - The linear power spectrum. 
Defaults to linear_power
    """

    def __init__(self, theta_min, theta_max, 
                 window_function_a, window_function_b, 
                 cosmo_dict=None, input_hod=None, halo_dict=None, 
                 powSpec=None, **kws):

        self.log_theta_min = numpy.log10(theta_min)
        self.log_theta_max = numpy.log10(theta_max)
        self.theta_array = numpy.logspace(
            self.log_theta_min, self.log_theta_max,
            defaults.default_precision["corr_npoints"])
        if theta_min==theta_max:
            self.log_theta_min = numpy.log10(theta_min)
            self.log_theta_max = numpy.log10(theta_min)
            self.theta_array = numpy.array([theta_min])
        self.wtheta_array= numpy.zeros(self.theta_array.size)

        # Hard coded, but we shouldn't expect halos outside of this range.
        self._k_min = 0.001
        self._k_max = 100.0

        if cosmo_dict == None:
            cosmo_dict = defaults.default_cosmo_dict
        self.cosmo_dict = cosmo_dict

        self.kernel = kernel.Kernel(self._k_min*theta_min, 
                                    self._k_max*theta_max,
                                    window_function_a, window_function_b,
                                    cosmo_dict)

        self.D_z = self.kernel.cosmo.growth_factor(self.kernel.z_bar)
                      
        if halo_dict is None:
            halo_dict = defaults.default_halo_dict
        self.halo = halo.HaloExclusion(input_hod, self.kernel.z_bar, 
                                       cosmo_dict, halo_dict)
        if powSpec==None:
            powSpec = 'linear_power'
        try:
            self.power_spec = self.halo.__getattribute__(powSpec)
        except AttributeError or TypeError:
            print "WARNING: Invalid input for power spectra varriable,"
            print "\t setting to linear_power"
            self.power_spec = self.halo.__getattribute__('linear_power')

    def set_redshift(self, redshift):
        self.kernel.z_bar = redshift
        self.D_z = self.kernel.cosmo.growth_factor(self.kernel.z_bar)
        self.halo.set_redshift(self.kernel.z_bar)
            
    def set_cosmology(self, cosmo_dict):
        self.kernel.set_cosmology(cosmo_dict)
        self.D_z = self.kernel.cosmo.growth_factor(self.kernel.z_bar)
        self.halo.set_cosmology(cosmo_dict, self.kernel.z_bar)

    def set_power_spectrum(self, powSpec):
        try:
            self.power_spec = self.halo.__getattribute__(powSpec)
        except AttributeError or TypeError:
            print "WARNING: Invalid input for power spectra varriable,"
            print "\t setting to 'linear_power'"
            self.power_spec = self.halo.__getattribute__('linear_power')

    def set_halo(self, halo_dict):
        self.halo.set_halo(halo_dict)

    def set_hod(self, input_hod):
        self.halo.set_hod(input_hod)

    def compute_correlation(self):
        for idx,theta in enumerate(self.theta_array):
            self.wtheta_array[idx] = self._correlation(theta)

    def _correlation(self, theta):
        ln_kmin = numpy.log(self._k_min)
        ln_kmax = numpy.log(self._k_max)
        wtheta = integrate.romberg(
            self._correlation_integrand, 
            ln_kmin, ln_kmax, args=(theta,), vec_func=True,
            tol=defaults.default_precision["corr_precision"])
        return wtheta

    def _correlation_integrand(self, ln_k, theta):
        return 1.0

    def write(self, output_file_name):
        f = open(output_file_name, "w")
        for theta, wtheta in zip(
            self.theta_array, self.wtheta_array):
            f.write("%1.10f %1.10f\n" % (theta/degToRad, wtheta))
        f.close()

class MagCorrelation(Correlation):
    def __init__(self, theta_min, theta_max, 
                 window_function_galaxy, window_function_convergence, 
                 cosmo_dict=None, input_hod=None, halo_dict=None, 
                 powSpec='power_gm', **kws):
        Correlation.__init__(self, theta_min, theta_max,
                             window_function_galaxy,
                             window_function_convergence,
                             cosmo_dict, input_hod, halo_dict, powSpec,
                             **kws)

    def _correlation_integrand(self, ln_k, theta):
        dln_k = 1.0
        k = numpy.exp(ln_k)
        dk = k*dln_k
        return 1.0/(2.0*numpy.pi)*(
            dk*k*self.power_spec(k)/(self.D_z*self.D_z)*
            self.kernel.kernel(numpy.log(k*theta)))

class AutoCorrelation(Correlation):

    def __init__(self, theta_min, theta_max, 
                 window_function_galaxy, 
                 cosmo_dict=None, input_hod=None, halo_dict=None, 
                 powSpec='power_gg', **kws):
        self.window_function = window_function_galaxy
        Correlation.__init__(self, theta_min, theta_max,
                             window_function_galaxy, window_function_galaxy,
                             cosmo_dict, input_hod, halo_dict, powSpec,
                             **kws)

    def _correlation_integrand(self, ln_k, theta):
        dln_k = 1.0
        k = numpy.exp(ln_k)
        dk = k*dln_k
        return 1.0/(2.0*numpy.pi)*(
            dk*k*self.power_spec(k)/(self.D_z*self.D_z)*
            self.kernel.kernel(numpy.log(k*theta)))

class GalaxyGalaxyLensing(Correlation):

    def __init__(self, theta_min, theta_max, 
                 window_function_galaxy, window_function_convergence, 
                 cosmo_dict=None, input_hod=None, halo_dict=None, 
                 powSpec='power_gm', **kws):
        self.log_theta_min = numpy.log10(theta_min)
        self.log_theta_max = numpy.log10(theta_max)
        self.theta_array = numpy.logspace(self.log_theta_min,
                                      self.log_theta_max, 50)
        if theta_min==theta_max:
            self.log_theta_min = numpy.log10(theta_min)
            self.log_theta_max = numpy.log10(theta_min)
            self.theta_array = numpy.array([theta_min])
        self.wtheta_array= numpy.zeros(self.theta_array.size)

        # Hard coded, but we shouldn't expect halos outside of this range.
        self._k_min = 0.001
        self._k_max = 100.0

        if cosmo_dict == None:
            cosmo_dict = defaults.default_cosmo_dict
        self.cosmo_dict = cosmo_dict

        self.kernel = kernel.GalaxyGalaxyLensingKernel(
            self._k_min*theta_min, 
            self._k_max*theta_max,
            window_function_galaxy, window_function_convergence,
            cosmo_dict)

        self.D_z = self.kernel.cosmo.growth_factor(self.kernel.z_bar)
                      
        if halo_dict is None:
            halo_dict = defaults.default_halo_dict
        self.halo = halo.HaloSatelliteCentral(input_hod, self.kernel.z_bar, 
                                              cosmo_dict, halo_dict)
        if powSpec==None:
            powSpec = 'linear_power'
        try:
            self.power_spec = self.halo.__getattribute__(powSpec)
        except AttributeError or TypeError:
            print "WARNING: Invalid input for power spectra varriable,"
            print "\t setting to linear_power"
            self.power_spec = self.halo.__getattribute__('linear_power')

    def _correlation_integrand(self, ln_k, theta):
        dln_k = 1.0
        k = numpy.exp(ln_k)
        dk = k*dln_k
        return 1.0/(4.0*numpy.pi)*(
            dk*k*self.power_spec(k)/(self.D_z*self.D_z)*
            self.kernel.kernel(numpy.log(k*theta)))
