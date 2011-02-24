""">>> param.py <<<

This module holds the classes for cosmological and other parameters,
provides default values  and an interface for accessing them.

Current revision:
    ID:         $Id: param.py 59 2008-02-28 20:36:04Z neyrinck $
    Date:       $Date: 2008-02-28 10:36:04 -1000 (Thu, 28 Feb 2008) $
    Revision:   $Revision: 59 $

(C) 2005 The CosmoPy Team (see Copyright for details)
"""

import numpy
import unittest
REVISION = "$Revision: 59 $"

# cosmological parameters have CAMB compatible names
defaultCosmoParamDict = {
        "ombh2": 0.02253,
        "omch2": 0.1122,
        "omnuh2": 0,
        "omk": 0,
        "hubble": 70.4,
        "w": -1,
        "cs2_lam": 1,
        "omega_baryon": 0.0455,
        "omega_cdm": 0.226,
        "omega_lambda": 0.728,
        "omega_neutrino": 0,
        "temp_cmb": 2.726,
        "helium_fraction": 0.249,
        "massless_neutrinos": 3.,
        "massive_neutrinos": 0,
        "re_redshift": 10.3,
        "re_ionization_frac": 1,
        "scalar_amp": [2.3e-9],
        "scalar_spectral_index": [0.967],
        "scalar_nrun": [0],
        "tensor_spectral_index": [0],
        "initial_ratio": [1]
        }

# note that all cosmological parameters are also Camb parameters
# initial_vector is a unique case, since it needs to print differently
# from other vector parameters, therefore I defined it is a string
# "1 0 0 0 0", but there surely must be a more elegant solution

defaultCambExtraParamDict = {
        "output_root": "tmpCosmoPyCamb",
        "get_scalar_cls": "F",
        "get_vector_cls": "F",
        "get_tensor_cls": "F",
        "get_transfer": "T",
        "do_lensing": "F",
        "do_nonlinear": 0,
        "l_max_scalar": 2000,
        "k_eta_max_scalar": 4000,
        "l_max_tensor": 1500,
        "k_eta_max_tensor": 3000,
        "reionization": "T",
        "re_use_optical_depth": "F",
        "initial_power_num": 1,
        "initial_condition": 1,
        "initial_vector": "-1 0 0 0 0",
        "vector_mode": 0,
        "COBE_normalize": "F",
        "CMB_outputscale": 7.4311e12,
        "transfer_high_precision": "F",
        "transfer_kmax": 100, # This should be set to ~100 for halo-model calcs,
                             # even if the intended range extends past 100
        "transfer_k_per_logint": 0,
        "transfer_num_redshifts": 1,
        "transfer_interp_matterpower": "T",
            # for now, use "T" for this new feature if you use pkInterp(k),
            # until non-uniform interpolation is supported.
        "transfer_redshift": [0],
        "transfer_filename": ["transfer_out.dat"],
        "transfer_matterpower": ["matterpower.dat"],
        "use_physical": "T",
        "scalar_output_file": "scalCls.dat",
        "vector_output_file": "vecCls.dat",
        "tensor_output_file": "tensCls.dat",
        "total_output_file": "totCls.dat",
        "lensed_output_file": "lensedCls.dat",
        "FITS_filename": "scalCls.fits",
        "feedback_level": 1,
        "lensing_method": 1,
        "accurate_BB": "F",
        "recombination": 1,
        "massive_nu_approx": 1,
        "accurate_polarization": "T",
        "accurate_reionization": "F",
        "do_tensor_neutrinos": "F",
        "do_late_rad_trunction": "T",
        "number_of_threads": 0,
        "accuracy_boost": 1,
        "l_accuracy_boost": 1,
        "l_sample_boost": 1
    }

defaultHaloModelParamDict = {
    "z":0.,  # redshift
    "delta":0.,  # Change to something else for, e.g. delta=1.686
    # Mass function parameters (param in Sheth Mo & Tormen, Eq.6)
    "st_big_a": 0.,  # PS: 0.5, ST: 0.322
    # To normalize mass function automatically (Int n(m,z) m dm = 1),
    # set st_big_a = 0.
    "stq": 0.3,  # PS: 0., ST:  0.3
    "st_little_a": 0.707,  # PS: 1., ST:  0.707 or 0.75 (SMT Eq. 5: sta = 1.)
    # Halo density profile parameters: parameterized in Cooray & Sheth
    "dpalpha": 1,  # NFW:1 Hernquist:1 Moore: 1.5 Moore(alphaform) = 1.5
    "dpbetanfw": 2,  # NFW:2 Hernquist:3 Moore: 1 Moore(alphaform) = 1.5
    "dpbetam": 1,  # NFW:1 Hernquist:1 Moore: 1.5 Moore(alphaform) = 1
    "nfw": "y",  # If "n", manually Fourier-transform halo profiles
                 # If "y", use Ci(x), Si(x) Numerical Recipes functions
                 # "analyticgamma": 0., gamma in analytic model.
                 # If zero, don't use analytic model.
                 # Standard in model: 1/6.  Deprecated.

    "k": 10.**numpy.arange(-2.,2.,0.1),
    "dofthp": 1,  # Calculate Fourier-transformed halo profiles (yes=1)?
                  # Optional since it can take a while, and doesn"t
                  # necessarily have to be redone.

    # concentration parameter distribution
    # default from Bullock et al. (2001) for NFW profile
    "cbarcoef": 9.,
    "cbarslope": -0.13,
    "sigmalogc": 0.2,  # Bullock et al: 0.2
    "concquadorder": 5,  # Order of Gaussian quadrature for concentration
                         # parameter

    # Implementation parameters
    "lomass":1e-156,  # in units of chimps (mass in Cubic h^-1 Mpcs)
                      # if 0, go to a scale based on the extent of c.pk
    "himass":3e8,
    "integrateto":0.,
    "massdivsperdex":10,
    "startki" : 100,    # Use this starting k index, relative to Camb
    "strideki": 10,     # the stride to take in k (camb samples really densely)
    "r_efoldings": 10,  # In density profile Fourier-transforming, the number
                        # of e-foldings of radius over which to integrate
    "numericalfthpaccuracy": 1e-6,  # log of ratio between values of r in d.p.
                                    # integration
    "rcutratio": 1.,    # ratio of halo cutoff radius to virial radius 
                        # (not yet implemented)
    "dcth": 0.1,        # interval in cos(theta) in angular averages
    "cth_endpts": "n",  # Allow endpoints in cos(theta) array?
    "dlna": 0.1,        # step in (natural) log of dlnP/dlnA derivative
    "largescalenormal": 0,  # whether (1) or not (0) to force the
                            # 2h term to be normalized to
                            # be unbiased wrt the linear power spectrum.
    "outputalltterms": 2,   # 0: just output total covmat in halocov
                            # 1: output ans,t1h,t2h,t3h,t4h
                            # 2: output ans,t1h, split other terms
    "use_sea_2h": 0,  # 0: use linear power spectrum for 2 halo term
                      # 1: use smith et al
    "whichp":"mm",    # do matter-matter statistics.
                      # "gg" (galaxy-galaxy) also acceptable.
                      # "gm" (galaxy-mass) also acceptable.

    # HOD parameters.  Notation as in Scoccimarro et al. (2001)
    "mcut_msunh": 4e12,  # cut-off mass for haloes to start having > 1 gal.
                         # Blue: 4e12 Msun/h.  Red: 2.5e12 Msuh/h
    "m11_msunh": 1e11,   # mass below which a halo has zero galaxies
                         # {0,    m < m11

    "m13_msunh": 1e13,
    "n0": 0.7,       # params in <N_gal|m> = {N0, m11 < m < mcut
    "alphaexp":0.8,  #                       {N0(m/mcut)^alphaexp, m > mcut

    "hodkrav":1,  # use Kravtsov et al."s (2004) HOD parameterization.
                  # 2: use possibly more accurate parameterization
                  # (not yet implemented)
    "k_mmin_msun": 1e11,
    "k_m1overmmin": 30.,
    "k_c": 0.045,
    "k_betas": 1.,
    }

class CosmoParams(dict):
    """base class for cosmological parameters"""
    def __init__(self,**pardict):
        super(CosmoParams,self).__init__(defaultCosmoParamDict)
        for k in pardict:
            if k in self: self[k] = pardict[k]
            else: raise cex.UnknownParameterError(k)
    def __getattr__(self,name):
        return self[name.lower()]
    def __str__(self):
        r = "parameters in CosmoParams:\n"
        for k in self:
            r += "\t-> %s = %s\n" % (k, str(self[k]))
        return r[:-1]


class CambParams(dict):
    """base class for camb (including cosmological) parameters"""
    def __init__(self,**pardict):
        super(CambParams,self).__init__(defaultCosmoParamDict)
        super(CambParams,self).__init__(defaultCambExtraParamDict)
        for k in pardict:
            if self.has_key(k):
                self[k] = pardict[k]
            else:
                raise cex.UnknownParameterError(k)
    def __getattr__(self,name):
        return self[name.lower()]
    def __str__(self):
        r = "parameters in CambParams:\n"
        for k in self:
            r += "\t-> %s = %s\n" % (k, str(self[k]))
        return r[:-1]


class HaloModelParams(dict):
    """base class for halo model parameters"""
    def __init__(self,**pardict):
        super(HaloModelParams,self).__init__(defaultHaloModelParamDict)
        for k in pardict:
            if k in self: self[k] = pardict[k]
            else: raise cex.UnknownParameterError(k)
    def __getattr__(self,name):
        return self[name.lower()]
    def __str__(self):
        r = "parameters in HaloModelParam:\n"
        for k in self:
            r += "\t-> %s = %s\n" % (k, str(self[k]))
        return r[:-1]

class UnknownParameterError(Exception):
    """exception raised when a parameter is not in the default dict"""
    pass

class CosmoParamsTests(unittest.TestCase):
    """unit tests for cosmological parameters"""
    def testParams0(self):
        self.cp = CosmoParams(hubble = 100)
        print self.cp
        self.failIf(self.cp.hubble != 100)

    def testParams1(self):
        try:
            self.cp = CosmoParams(huble = 100)
            self.failIf(True)
        except UnknownParameterError:
            self.failIf(False)

class CambParamsTests(unittest.TestCase):
    """unit tests for cosmological parameters"""
    def testParams0(self):
        self.cp = CambParams(hubble = 100,output_root="xxx")
        print self.cp
        self.failIf(self.cp.hubble != 100 or self.cp.output_root != "xxx")

    def testParams1(self):
        try:
            self.cp = CosmoParams(outputroot = "xxx")
            self.failIf(True)
        except UnknownParameterError:
            self.failIf(False)


def main():
    unittest.main()

if __name__=="__main__":
    main()
