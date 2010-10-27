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

__author__ = "Chris Morrison <morrison.chrisb@gmail.com>"

class Correlation(object):

    """
    Template class for correlation function, theta should be in radians. The integrand is over ln_k_theta.
    """

    def __init__(self, theta_min, theta_max, window_function_a, window_function_b, camb_param=None, input_hod=None, halo_param=None, powSpec=None, **kws):

        self.log_theta_min = np.log10(theta_min)
        self.log_theta_max = np.log10(theta_max)
        self.theta_array = np.logspace(self.log_theta_min, 
                                       self.log_theta_max, 20)
        self.wtheta_array= np.zeros(self.theta_array.size)

        # Hard coded, but we shouldn't expect halos outside of this range.
        self._k_min = 0.001
        self._k_max = 100.0

        self.kernel = kernel.Kernel(self._k_min*theta_min, self._k_max*theta_max,
                                    window_function_a, window_function_b)
        if camb_param == None:
            camb_param = param.CambParams(**kws)
        self.kernel.set_cosmology(camb_param)
                      
        if halo_param is None:
            halo_param = param.HaloModelParams(**kws)
        self.halo = halo.Halo(input_hod, self.kernel.z_bar, camb_param, halo_param)
            
    def set_cosmology(self, camb_param):
        self.kernel.set_cosmology(camb_param)
        self.halo.set_cosmology(camb_param, self.kernel.z_bar)

    def set_halo(self, halo_param):
        self.halo.set_halo(halo_param)

    def set_hod(self, input_hod):
        self.halo.set_hod(input_hod)

    def compute_correlation(self):
        for idx in xrange(self.theta_array.size):
            print "Computing Theta at", self.theta_array[idx]
            wtheta, wtheta_err = integrate.quad(self._correlation_integrand, 
                                                self.kernel.ln_ktheta_min,
                                                self.kernel.ln_ktheta_max,
                                                args=(self.theta_array[idx]),
                                                limit=200)
            self.wtheta_array[idx] = wtheta

    def _correlation_integrand(self, ln_ktheta, theta):
        k = np.exp(ln_ktheta)/theta
        dk = k
        return dk*k*self.halo.power_gg(k)*self.kernel.kernel(ln_ktheta)

    def write(self, output_file_name):
        f = open(output_file_name, "w")
        for theta, wtheta in zip(
            self.theta_array, self.wtheta_array):
            f.write("%1.10f %1.10f\n" % (theta, wtheta))
        f.close()
