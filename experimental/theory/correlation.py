import camb
import cosmology
import halo
import hod
import kernel
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
                 camb_param=None, input_hod=None, halo_param=None, 
                 powSpec=None, **kws):

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

        self.kernel = kernel.Kernel(self._k_min*theta_min, 
                                    self._k_max*theta_max,
                                    window_function_a, window_function_b)
        if camb_param == None:
            camb_param = param.CambParams(**kws)
        self.kernel.set_cosmology(camb_param)

        self.D_z = self.kernel.cosmo.growth_factor(self.kernel.z_bar)
                      
        if halo_param is None:
            halo_param = param.HaloModelParams(**kws)
        self.halo = halo.Halo(input_hod, self.kernel.z_bar, camb_param, 
                              halo_param)
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
            
    def set_cosmology(self, camb_param):
        self.kernel.set_cosmology(camb_param)
        self.D_z = self.kernel.cosmo.growth_factor(self.kernel.z_bar)
        self.halo.set_cosmology(camb_param, self.kernel.z_bar)

    def set_power_spectrum(self, powSpec):
        try:
            self.power_spec = self.halo.__getattribute__(powSpec)
        except AttributeError or TypeError:
            print "WARNING: Invalid input for power spectra varriable,"
            print "\t setting to 'linear_power'"
            self.power_spec = self.halo.__getattribute__('linear_power')

    def set_halo(self, halo_param):
        self.halo.set_halo(halo_param)

    def set_hod(self, input_hod):
        self.halo.set_hod(input_hod)

    def compute_correlation(self):
        for idx,theta in enumerate(self.theta_array):
            self.wtheta_array[idx] = self._correlation(theta)

    def _correlation(self, theta):
        ln_kmin = numpy.log(self._k_min)
        ln_kmax = numpy.log(self._k_max)
        wtheta, wtheta_err = integrate.quad(self._correlation_integrand, 
                                            ln_kmin,
                                            ln_kmax,
                                            args=(theta,),
                                            limit=200)
        print "Theta:",theta,"Wtheta:",wtheta
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
                 window_function_a, window_function_b, 
                 camb_param=None, input_hod=None, halo_param=None, 
                 powSpec='power_gm', **kws):
        Correlation.__init__(self, theta_min, theta_max,
                             window_function_a, window_function_b,
                             camb_param, input_hod, halo_param, powSpec,
                             **kws)

    def _correlation_integrand(self, ln_k, theta):
        dln_k = 1.0
        k = numpy.exp(ln_k)
        dk = k*dln_k
        return 2.0*numpy.sqrt(2)/(8.0*numpy.pi)*(
            dk*k*self.power_spec(k)/(self.D_z*self.D_z)*
            self.kernel.kernel(numpy.log(k*theta)))

class AutoCorrelation(Correlation):

    def __init__(self, theta_min, theta_max, 
                 window_function, 
                 camb_param=None, input_hod=None, halo_param=None, 
                 powSpec='power_gg', **kws):
        self.window_function = window_function
        Correlation.__init__(self, theta_min, theta_max,
                             window_function, window_function,
                             camb_param, input_hod, halo_param, powSpec,
                             **kws)

    def _correlation_integrand(self, ln_k, theta):
        dln_k = 1.0
        k = numpy.exp(ln_k)
        dk = k*dln_k
        return 2.0*numpy.sqrt(2)/(8.0*numpy.pi)*(
            dk*k*self.power_spec(k)/(self.D_z*self.D_z)*
            self.kernel.kernel(numpy.log(k*theta)))

class AutoCorrelationDeltaFunction(Correlation):

    def __init__(self, theta_min, theta_max, redshift, input_cosmology=None,
               input_hod=None, input_halo=None, powSpec=None, **kws):
        self.log_theta_min = numpy.log10(theta_min)
        self.log_theta_max = numpy.log10(theta_max)
        self.theta_array = numpy.logspace(self.log_theta_min, 
                                       self.log_theta_max, 20)
        self.wtheta_array= numpy.zeros(self.theta_array.size)
        self.z = redshift
        self._k_min = 0.001
        self._k_max = 100.0
        
        self.ln_ktheta_min = numpy.log(self._k_min*theta_min)
        self.ln_ktheta_max = numpy.log(self._k_max*theta_max)

        dln_ktheta = (self.ln_ktheta_max - self.ln_ktheta_min)/200.0
        self._ln_ktheta_array = numpy.arange(
            self.ln_ktheta_min, self.ln_ktheta_max + dln_ktheta, dln_ktheta)

        if input_cosmology == None:
            self.cosmology = cosmology.SingleEpoch(redshift,
                                                   param.CambParams(**kws))
        else:
            self.cosmology = input_cosmology
        self.halo = input_halo
            

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
        return (dk*k*self.power_spec(k)*
                special.j0(k*theta*cosmo.comoving_distance()))

class CorrelationDeltaFunction(Correlation):
    
    """
    Correlation function for delta functions assumed for the redshift 
    distribution of both populations being correlated.
    """

    def __init__(self, theta_min, theta_max, redshift_a, redshift_b, 
                 camb_param=None, input_hod=None, halo_param=None, 
                 powSpec=None, **kws):
        self.log_theta_min = numpy.log10(theta_min)
        self.log_theta_max = numpy.log10(theta_max)
        self.theta_array = numpy.logspace(self.log_theta_min, 
                                       self.log_theta_max, 20)
        if theta_min == theta_max:
            self.log_theta_min = numpy.log10(theta_min)
            self.log_theta_max = numpy.log10(theta_min)
            self.theta_array = numpy.array([theta_min])
        self.wtheta_array= numpy.zeros(self.theta_array.size)

        self.z_a = redshift_a
        self.z_b = redshift_b

        # Hard coded, but we shouldn't expect halos outside of this range.
        self._k_min = 0.001
        self._k_max = 100.0
        
        self.ln_ktheta_min = numpy.log(self._k_min*theta_min)
        self.ln_ktheta_max = numpy.log(self._k_max*theta_max)
        dln_ktheta = (self.ln_ktheta_max - self.ln_ktheta_min)/200.0
        self._ln_ktheta_array = numpy.arange(
            self.ln_ktheta_min, self.ln_ktheta_max + dln_ktheta, dln_ktheta)

        if camb_param == None:
            camb_param = param.CambParams(**kws)
        self.cosmology = cosmology.MultiEpoch(0,self.z_b+self.z_a,camb_param)
        self.xiA_a = self.cosmology.comoving_distance(self.z_a)
        self.xiA_b = self.cosmology.comoving_distance(self.z_b)
        print (self.xiA_a, self.xiA_b, self.cosmology.h, 
               self.cosmology.omega_m0, 1+self.z_a, self.xiA_b-self.xiA_a,
               self.xiA_b*self.xiA_a,
               self.cosmology.growth_factor(self.z_a))
        self.cosmology_prefactor = (3.0/2.0*(self.cosmology.h*100/
                                             speed_of_light)**2*
                                    self.cosmology.omega_m0*
                                    (1+self.z_a)*
                                    (self.xiA_b-self.xiA_a)/
                                    (self.xiA_b*self.xiA_a))
        print self.cosmology_prefactor

        self.ln_s_min = numpy.log(self._k_min*self.xiA_a)
        self.ln_s_max = numpy.log(self._k_max*self.xiA_b)
        dln_s = (self.ln_s_max - self.ln_s_min)/200.0
        self._ln_s_array = numpy.arange(
            self.ln_s_min, self.ln_s_max + dln_s, dln_s)
                      
        if halo_param is None:
            halo_param = param.HaloModelParams(**kws)
        self.halo = halo.Halo(input_hod, self.z_a, camb_param, halo_param)

    def compute_correlation(self):
        for idx in xrange(self.theta_array.size):
            print "Computing Theta at", self.theta_array[idx]
            wtheta, wtheta_err = integrate.quad(self._correlation_integrand, 
                                                self.ln_s_min,
                                                self.ln_s_max,
                                                args=(self.theta_array[idx]),
                                                limit=200)
            self.wtheta_array[idx] = self.cosmology_prefactor*wtheta

    def _correlation_integrand(self, ln_s, theta):
        dln_s = 1.0
        s = numpy.exp(ln_s)
        ds = s*dln_s
        k = s/self.xiA_a
        return (ds*s*self.halo.power_mm(k)
                *special.j0(s*theta))

    def write_kernel(self, output_file_name):
        f = open(output_file_name, "w")
        theta = 1.0/60.0*numpy.pi/180.0
        for ln_s in self._ln_s_array:
            s = numpy.exp(ln_s)
            k = s/self.xiA_a
            f.write("%1.10f %1.10f %1.10f\n" % (k, self._correlation_integrand(
                        ln_s, theta)/s**2,special.j0(s*theta)))
        f.close()

