from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import integrate as In
from scipy import special as S
import param  # cosmological parameter object from Cosmopy
import cosmology
import numpy
import copy
import cosmology

"""This is a set of classes for constructing an angular correlation kernel.

To calculate the Limber's approximation of an angular correlation function,
we need a kernel object that integrates over our window functions and translates
between physical and angular scale.  These window functions can be simply the
redshift distribution of our galaxies, a lensing function, an ISW potential
decay function and so on.  The idea here is that our kernel is a generic object
that takes two window function objects and performs all of the necessary
integrals.  It can also return the peak in the redshift sensitivity, so we can
make the best approximation as to the appropriate power spectrum.
"""

__author__ = "Ryan Scranton <ryan.scranton@gmail.com>"

class dNdz(object):
    """Base class for a simple redshift distribution.

    This class handles all of the details of normalization and interpolation.
    Derived classes should be used for specific redshift distributions.
    """
    def __init__(self, z_min, z_max):
        self.z_min = z_min
        self.z_max = z_max
        dz = (self.z_max - self.z_min)/100
        self.norm = 1.0

    def normalize(self):
        norm, norm_err = In.quad(self.dndz, self.z_min, self.z_max)
        print "dN/dz normalization = %1.10f, err = %1.10f" % (norm, norm_err)

        self.norm = 1.0/norm

    def raw_dndz(self, redshift):
        return 1.0

    def dndz(self, redshift):
        if redshift <= self.z_max and redshift >= self.z_min:
            return self.norm*self.raw_dndz(redshift)
        else:
            return 0.0

class dNdzGaussian(dNdz):
    """Derived class for a Gaussian-shaped redshift distribution.

    dNdz ~ exp(-(z-z0)^2/sigma_z^2)

    """
    def __init__(self, z_min, z_max, z0, sigma_z):
        dNdz.__init__(self, z_min, z_max)
        self.z0 = z0
        self.sigma_z = sigma_z

    def raw_dndz(self, redshift):
        return numpy.exp(-1.0*(redshift-self.z0)*(redshift-self.z0)/
                         (2.0*self.sigma_z*self.sigma_z))

class dNdzMagLim(dNdz):
    """Derived class for a magnitude-limited redshift distribution.

    dNdz ~ z^a*exp(-(z/z0)^b)

    """
    def __init__(self, z_min, z_max, a, z0, b):
        dNdz.__init__(self, z_min, z_max)
        self.a = a
        self.z0 = z0
        self.b = b

    def raw_dndz(self, redshift):
        return ((redshift**self.a)*numpy.exp(-1.0*redshift/self.z0**self.b))

class WindowFunction(object):
    """Base class for an angular correlation window function.

    This object represents the window function for one of the two fields going
    into a correlation measurement, expressed as a function of comoving
    distance.  The details of the window function depends on the field
    involved (galaxy distribution, lensing potential, ISW potential, etc.), but
    the base class defines the API which is necessary in order to integrate the
    window function over comoving distance.

    In general, a proper calculation of the window function may be expensive,
    involving integrals of its own, so we want to evaluate the window function
    at a set number of points and then spline over them when it comes time to
    integrate over the wave function.  The former is done via the
    raw_window_function() method, which will be re-implemented by each of the
    derived classes and then the kernel object will sample the spline via the
    window_function() method.
    """
    def __init__(self, z_min, z_max, camb_param=None, **kws):
        self.initialized_spline = False

        self.z_min = z_min
        self.z_max = z_max

        if camb_param is None:
            camb_param = param.CambParams(**kws)
        self.set_cosmology(camb_param)

        self.chi_min = self.cosmo.comoving_distance(z_min)
        self.chi_max = self.cosmo.comoving_distance(z_max)

        dchi = (self.chi_max - self.chi_min)/200.0
        self._chi_array = numpy.arange(self.chi_min, self.chi_max + dchi, dchi)
        self._wf_array = numpy.zeros_like(self._chi_array)

    def set_cosmology(self, camb_param):
        self.cosmo = cosmology.MultiEpoch(self.z_min, self.z_max, camb_param)

        if self.initialized_spline:
            self._initialize_spline()

    def _initialize_spline(self):
        for idx in xrange(self._chi_array.size):
            self._wf_array[idx] = self.raw_window_function(self._chi_array[idx])
        self._wf_spline = InterpolatedUnivariateSpline(self._chi_array,
                                                       self._wf_array)
        self.initialized_spline = True

    def raw_window_function(self, chi):
        return 1.0

    def window_function(self, chi):
        if not self.initialized_spline:
            self._initialize_spline()

        if chi <= self.chi_max and chi >= self.chi_min:
            return self._wf_spline(chi)[0]
        else:
            return 0.0

    def write(self, output_file_name):
        f = open(output_file_name, "w")
        for chi, wf in zip(self._chi_array, self._wf_array):
            f.write("%1.10f %1.10f\n" % (chi, wf))
        f.close()

class WindowFunctionGalaxy(WindowFunction):
    """WindowFunction class for a galaxy distribution.

    This derived class takes the standard WindowFunction arguments along with
    a redshift distribution and turns it into a proper WindowFunction for
    kernel integration:

    W(chi) = dN/dz dz/dchi

    for comoving distance chi.
    """
    def __init__(self, redshift_dist, camb_param=None, **kws):
        self._redshift_dist = redshift_dist
        self._redshift_dist.normalize()

        WindowFunction.__init__(self, redshift_dist.z_min, redshift_dist.z_max,
                                camb_param, **kws)

    def raw_window_function(self, chi):
        z = self.cosmo.redshift(chi)

        dzdchi = 1.0/self.cosmo.E(z)

        return dzdchi*self._redshift_dist.dndz(z)

class WindowFunctionConvergence(WindowFunction):
    """WindowFunction class for magnification of a background sample.

    This derived class calculates the magnification effect of a background
    sample as a function of comoving distance chi.  In essence, given a sample
    with redshift distribution dN/dz, what is the weighted fraction of that
    sample that is beyond chi:

    g(chi) = chi*int(chi, inf, dN/dz dz/dchi' (1.0 - chi/chi'))

    and the window function is

    W(chi) = 3*omega_m*g(chi)/a

    where we have omitted the (2.5*s - 1) factor that is due to the number count
    slope of the sample.
    """
    def __init__(self, redshift_dist, camb_param=None, **kws):
        self._redshift_dist = redshift_dist
        self._redshift_dist.normalize()

        self._g_chi_min = 0.0
        # Even though the input distribution may only extend between some bounds
        # in redshift, the lensing kernel will extend across z = [0, z_max)
        WindowFunction.__init__(self, 0.0, redshift_dist.z_max,
                                camb_param, **kws)
        self._g_chi_min = (
            self.cosmo.comoving_distance(self._redshift_dist.z_min))

    def raw_window_function(self, chi):
        a = 1.0/(1.0 + self.cosmo.redshift(chi))

        chi_bound = chi
        if chi_bound < self._g_chi_min: chi_bound = self._g_chi_min

        g_chi, g_chi_err = In.quad(self._lensing_integrand, chi_bound,
                                   self.chi_max, args=(chi,))

        g_chi *= self.cosmo.H0*self.cosmo.H0*chi

        return 3.0*self.cosmo.omega_m0*g_chi/a

    def _lensing_integrand(self, chi, chi0):
        z = self.cosmo.redshift(chi)

        dzdchi = 1.0/self.cosmo.E(z)

        return dzdchi*self._redshift_dist.dndz(z)*(chi - chi0)/chi

class Kernel(object):
    """Container class for calculating correlation function kernels.

    A kernel is an integtral over the product of two window functions
    representing the spatial extent of two fields (or one field for the case of
    an autocorrelation).  In addition, there is a Bessel function which
    incorporates the projected angular dependence of the correlation function.
    This means that the kernel is a function of k*theta, where k is in h/Mpc
    and theta is in radians:

    K(k, theta) = 4pi^2*int(0, inf, D^2(chi)*W_a(chi)*W_b(chi)*J_0(k*theta*chi))

    In addition to providing the kernel function, a kernel object also
    calculates z_bar, the peak in the kernel redshift sensitivity.
    """
    def __init__(self, ktheta_min, ktheta_max,
                 window_function_a, window_function_b,
                 camb_param=None, **kws):
        self.initialized_spline = False

        self.ln_ktheta_min = numpy.log(ktheta_min)
        self.ln_ktheta_max = numpy.log(ktheta_max)

        if camb_param is None:
            camb_param = param.CambParams(**kws)

        self.window_function_a = window_function_a
        self.window_function_b = window_function_b

        self.window_function_a.set_cosmology(camb_param)
        self.window_function_b.set_cosmology(camb_param)

        self.chi_min = self.window_function_a.chi_min
        self.z_min = self.window_function_a.z_min
        if self.window_function_b.chi_min < self.chi_min:
            self.chi_min = self.window_function_b.chi_min
            self.z_min = self.window_function_b.z_min

        self.chi_max = self.window_function_a.chi_max
        self.z_max = self.window_function_a.z_max
        if self.window_function_b.chi_max > self.chi_max:
            self.chi_max = self.window_function_b.chi_max
            self.z_max = self.window_function_b.z_max
        self.cosmo = cosmology.MultiEpoch(self.z_min, self.z_max, camb_param)

        dln_ktheta = (self.ln_ktheta_max - self.ln_ktheta_min)/200.0
        self._ln_ktheta_array = numpy.arange(
            self.ln_ktheta_min, self.ln_ktheta_max + dln_ktheta, dln_ktheta)
        self._kernel_array = numpy.zeros_like(self._ln_ktheta_array)

        self._find_z_bar()

    def _find_z_bar(self):
        dz = (self.z_max - self.z_min)/100

        kernel_max = -1.0e30
        for z in numpy.arange(self.z_min, self.z_max, dz):
            kernel = self._kernel_integrand(
                self.cosmo.comoving_distance(z), 1.0)
            if kernel > kernel_max:
                self.z_bar = z

    def _initialize_spline(self):
        kernel_max = -1.0e30
        calculate_kernel = True
        ratio_limit = 1.0e-4
        for idx in xrange(self._ln_ktheta_array.size):
            kernel = 0.0
            if calculate_kernel:
                kernel = self.raw_kernel(self._ln_ktheta_array[idx])

            if numpy.abs(kernel) > kernel_max:
                kernel_max = numpy.abs(kernel)
            if numpy.abs(kernel/kernel_max) < ratio_limit:
                calculate_kernel = False
            self._kernel_array[idx] = kernel
            print "%1.10f %1.10f" % (self._ln_ktheta_array[idx],
                                     self._kernel_array[idx])

        self._kernel_spline = InterpolatedUnivariateSpline(
            self._ln_ktheta_array, self._kernel_array)

        self.initialized_spline = True

    def set_cosmology(self, camb_param):
        self.initialized_spline = False

        self.window_function_a.set_cosmology(camb_param)
        self.window_function_b.set_cosmology(camb_param)
        
        self.chi_min = self.window_function_a.chi_min
        self.z_min = self.window_function_a.z_min
        if self.window_function_b.chi_min < self.chi_min:
            self.chi_min = self.window_function_b.chi_min
            self.z_min = self.window_function_b.z_min

        self.chi_max = self.window_function_a.chi_max
        self.z_max = self.window_function_a.z_max
        if self.window_function_b.chi_max > self.chi_max:
            self.chi_max = self.window_function_b.chi_max
            self.z_max = self.window_function_b.z_max
        self.cosmo = cosmology.MultiEpoch(self.z_min, self.z_max, camb_param)

        self._find_z_bar()  

    def raw_kernel(self, ln_ktheta):
        ktheta = numpy.exp(ln_ktheta)

        kernel, kernel_err = In.quad(self._kernel_integrand, 1.01*self.chi_min,
                                     0.99*self.chi_max, args=(ktheta,),
                                     limit=200)

        return 4.0*numpy.pi*numpy.pi*kernel

    def _kernel_integrand(self, chi, ktheta):
        D_z = self.window_function_a.cosmo.growth_factor(
            self.window_function_a.cosmo.redshift(chi))

        return (self.window_function_a.window_function(chi)*
                self.window_function_b.window_function(chi)*
                D_z*D_z*S.j0(ktheta*chi))

    def kernel(self, ln_ktheta):
        if not self.initialized_spline:
            self._initialize_spline()

        if (ln_ktheta <= self.ln_ktheta_max and
            ln_ktheta >= self.ln_ktheta_min):
            return self._kernel_spline(ln_ktheta)[0]
        else:
            return 0.0

    def write(self, output_file_name):
        if not self.initialized_spline:
            self._initialize_spline()

        f = open(output_file_name, "w")
        for ln_ktheta, kernel in zip(
            self._ln_ktheta_array, self._kernel_array):
            f.write("%1.10f %1.10f\n" % (numpy.exp(ln_ktheta), kernel))
        f.close()
