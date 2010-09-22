from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import integrate
import param  # cosmological parameter object from Cosmopy
import numpy
import cosmology
import os
import types

"""Wrapper class for interfacing with the CAMB code.

Given a set of parameters, we want to use the CAMB code to generate a linear
power spectrum we can output for arbitrary wavenumber.
"""

__author__ = "Ryan Scranton <ryan.scranton@gmail.com"

class CambWrapper(object):
    """Wrapper for calling the CAMB code."""
    def __init__(self, camb_param=None, halo_param=None, **kws):
        self.camb_path = ""
        self.parameter_file = "camb_param.ini"

        if camb_param is None:
            camb_param = param.CambParams(**kws)

        self._camb_param = camb_param

        self._initialized = False

    def set_redshift(self, redshift):
        self._camb_param.transfer_redshift[0] = redshift

    def redshift(self):
        return self._camb_param.transfer_redshift[0]

    def run(self, parameter_file=None, cleanup_files=True):
        self._initialized = False

        camb_executable = os.path.join(self.camb_path, "camb")
        if not os.path.exists(camb_executable):
            print "%s does not exist!" % camb_executable
            return False

        if not parameter_file is None:
            self.parameter_file = parameter_file
        self.write_parameter_file(self.parameter_file)

        os.system(camb_executable + " " + self.parameter_file)

        power_spectrum_file = self.camb_param.output_root+"_matterpower.dat"
        if not os.path.exists(power_spectrum_file):
            print "%s does not exist!" % power_spectrum_file

        k_array, power_array = self.read_power_spectrum(power_spectrum_file)
        self.k_min = 1.01*k_array[0]
        self.k_max = 0.99*k_array[-1]

        self._log_k_array = numpy.zeros_like(k_array)
        self._log_power_array = numpy.zeros_like(k_array)

        for idx in xrange(k_array.size):
            self._ln_k_array[idx] = numpy.log(k_array[idx])
            self._ln_power_array[idx] = numpy.log(power_array[idx])
        self._ln_power_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, self._ln_power_array)

        if clean_files:
            self.cleanup_files()

        self._initialized = True

        return self._initialized

    def write_parameter_file(self, parameter_file):
        output_file = file(parameter_file, "w")

        for key in self._camb_param.iterkeys():
            if type(self._camb_param[key]) == ListType:
                for i,val in enumerate(self._camb_param[key]):
                    output_file.write(key+"("+str(i+1)+") = "+str(val)+"\n")
            else:
                output_file.write(key+" = "+str(self._camb_param[key])+"\n")

        output_file.close()

    def read_power_spectrum(self, power_spectrum_file):
        power_matrix = numpy.loadtxt(power_spectrum_file)

        return power_matrix[:,0], power_matrix[:,1]

    def cleanup_files(self):
        if os.path.exists(self.paramter_file):
            os.system("rm %s" % parameter_file)

        power_spectrum_file = self.camb_param.output_root+"_matterpower.dat"
        if os.path.exists(power_spectrum_file):
            os.system("rm %s" % power_spectrum_file)

        transfer_file = self.camb_param.output_root+"_transfer_out.dat"
        if os.path.exists(transfer_file):
            os.system("rm %s" % transfer_file)

        cl_file = self.camb_param.output_root+"_scalCls.dat"
        if os.path.exists(cl_file):
            os.system("rm %s" % cl_file)

    def linear_power(self, k):
        if not self._initialized:
            self.run()

        if k < self.k_max and k > self.k_min:
            return numpy.exp(self._ln_power_spline(numpy.log(k)))
        else:
            return 0.0

    def sigma_r(self, scale):
        sigma2, sigma2_err = integrate.quad(
            self._sigma_integrand, numpy.log(self.k_min),
            numpy.log(self.k_max), args=(scale,))
        sigma2 *= numpy.log(10.0)/(2.0*numpy.pi**2)

        return numpy.sqrt(sigma2)

    def _sigma_integrand(self, logk, scale):
        k = 10**logk
        kR = scale*k

        W = 9.0*(numpy.sin(kR) - kR*numpy.cos(kR))**2/(kR**6)

        return self.linear_power(k)*W*k**3

    def normalize(self, sigma8):
        current_sigma8 = self.sigma_r(8.0)

        for idx in xrange(self._ln_power_array):
            self._ln_power_array[idx] *= sigma8/current_sigma8
        self._ln_power_spline = InterpolatedUnivariateSpline(
            self._ln_k_array, self._ln_power_array)


