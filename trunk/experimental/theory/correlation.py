import camb
import cosmology
import halo
import hod
import kernel
import param
import numpy as np
from scipy import integrate

"""Classes describing different correlation functions

Given a kernel function and halo model the class should produce a correlation as function of theta. The correlation should specifily the range of theta to be computed.

"""

class Correlation(object):

    """
    Template class for correlation function, theta should be in radians. The integrand is over ln_k_theta.
    """

    def __init__(self, theta_min, theta_max, window_function_a, window_function_b, camb_param=None, hod_param=None, **kws):
        self.log_theta_min = np.log10(theta_min)
        self.log_theta_max = np.log10(theta_max)
        self.theta_array = np.logspace(self.log_theta_min, 
                                       self.log_theta_max, 
                                       20)
        self.wtheta_array= np.zeros(self.log_theta_array.size)

        self.window_function_a = window_function_a
        self.window_function_b = window_function_b

        self.z_min = self.window_function_a.z_min
        if self.z_min > self.window_function_b.z_min:
            self.z_min = self.window_function_b.z_min
        self.z_max = self.window_function_a.z_max
        if self.z_max > self.window_function_b.z_max:
            self.z_max = self.window_function_b.z_max
        
        if camb_param is None:
            camb_param = param.CambParams(**kws)
        self.cosmology = cosmology.MultiEpoch(self.z_min, self.z_max, camb_param)
        if hod_param is None:
            hod_param = param.HodModelParams(**kws)
        self.hod = hod.HODKravtsov(hod_param)
            
        

    def compute_correlation(self):
        for idx in xrange(self.theta_array.size):
            wtheta, wtheta_err = integrante.quad(self._correlation_integrand, 
                                                 self.kernel.ln_ktheta_min,
                                                 self.kernel.ln_ktheta_max,
                                                 args=(self.theta_array[idx]])
            self.wtheta_array[idx] = wtheta

    def _correlation_integrand(self, ln_k_theta, log_theta):
        k = np.exp(ln_k_theta)/theta
        dk = k
        return 

    def write(self, output_file_name):
        f = open(output_file_name, "w")
        for theta, wtheta in zip(
            self.theta_array, self.wtheta_array):
            f.write("%1.10f %1.10f\n" % (theta, wtheta))
        f.close()
