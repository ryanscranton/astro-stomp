from matplotlib import pyplot
import numpy

degToRad = numpy.pi/180.0

def cosmology_unit_test():
    import cosmology
    print "/n****************************"
    print "*                          *"
    print "* Testing Cosmology Module *"
    print "*                          *"
    print "****************************\n\n"
    print "Testing Single Epoch"
    print "****************************"

    ### Create single epoch cosmologies at the redshift specified
    ### outputs to stdout
    print "\tRedshift z=0.0"
    cosmo = cosmology.SingleEpoch(redshift=0.0)
    cosmo.write()

    print "\tRedshift z=0.5"
    cosmo.set_redshift(redshift=0.5)
    cosmo.write()

    print "\tRedshift z=1.0"
    cosmo.set_redshift(redshift=1.0)
    cosmo.write()
    
    print "\tRedshift z=2.0"
    cosmo.set_redshift(redshift=2.0)
    cosmo.write()

    print "\tRedshift z=3.0"
    cosmo.set_redshift(redshift=3.0)
    cosmo.write()

    ### Compute example multi epoch cosmologies from redshift z=0.0 to z=5.0
    ### plotted are the computed distnaces as a funciton of redshift and
    ### several other cosmological variables (Omega_m(z), Omega_L(z), etc.)
    print "\nTesting Multi Epoch"
    print "****************************"
    cosmo = cosmology.MultiEpoch(z_min=0.0, z_max=5.0)
    pyplot.figure(figsize=(10,8))
    pyplot.semilogy(cosmo.comoving_distance(cosmo._z_array),label=r"$D_M$")
    pyplot.semilogy(cosmo.angular_diameter_distance(cosmo._z_array),
                    label=r"$D_A$")
    pyplot.semilogy(cosmo.luminosity_distance(cosmo._z_array),
                    label=r"$D_L$")
    pyplot.legend(loc=0)
    pyplot.xlabel('z')
    pyplot.ylabel('Mpc/h')
    pyplot.savefig('test_cosmology_distance.pdf')
    
    cosmo = cosmology.MultiEpoch(0.0, 5.0)
    pyplot.figure(figsize=(10,8))
    pyplot.semilogy(cosmo._z_array, 
                    cosmo.growth_factor(cosmo._z_array),label=r"Growth(z)")
    pyplot.semilogy(cosmo._z_array,
                    cosmo.omega_m(cosmo._z_array),
                    label=r"$\Omega_m$")
    pyplot.semilogy(cosmo._z_array,
                    cosmo.omega_l(cosmo._z_array),
                    label=r"$\Omega_{\Lambda}$")
    pyplot.semilogy(cosmo._z_array, 
                    cosmo.delta_c(cosmo._z_array),
                    label=r"$\delta_c(z)$")
    pyplot.semilogy(cosmo._z_array, 
                    cosmo.sigma_r(8.0, cosmo._z_array),
                    label=r"$\sigma_8(z)$")
    pyplot.legend(loc=0)
    pyplot.xlabel('z')
    pyplot.savefig('test_cosmology_misc.pdf')
    

def mass_function_unit_test():
    import mass_function
    print "/n********************************"
    print "*                              *"
    print "* Testing Mass Function Module *"
    print "*                              *"
    print "********************************\n\n"
    
    ### Compute the mass function at redshift z=0.0
    mass_func = mass_function.MassFunction(redshift=0.0)
    mass_array = numpy.exp(mass_func._ln_mass_array)

    ### Plot the normalized overdensity nu, the mass function f_m, and the
    ### halo bias.
    pyplot.figure(figsize=(10, 8))
    pyplot.loglog(mass_array, mass_func.nu(mass_array), label=r"$\nu(M)$")
    pyplot.loglog(mass_array, mass_func.f_m(mass_array), label="dN/dM")
    pyplot.loglog(mass_array, mass_func.bias_m(mass_array), label="b(M)")
    pyplot.legend(loc=0)
    pyplot.xlabel('Mass [$M_{\odot}/h^2$]')
    pyplot.savefig('test_mass_function.pdf')

def halo_unit_test():
    import halo
    import hod
    print "/n***********************"
    print "*                     *"
    print "* Testing Halo Module *"
    print "*                     *"
    print "***********************\n\n"
    

    ### We test each of the 4 avalible power spectra as a function of redshift
    ### first create a halo occupation distribution object using a Zheng07 HOD
    zheng = hod.HODZheng(10**13.0, 0.15, 10**13.0, 10**14.0, 1.0)
    ### initialize the halo model object at z=0.0 with the Zheng HOD
    h = halo.Halo(redshift=0.0, input_hod=zheng)
    ### define the k space range to plot
    k_array = numpy.logspace(-3, 2, 200)
    pyplot.figure(figsize=(10,8))
    ### We all power spectrum as the unitless Delta^2(k)-P(k)*k**3
    ### Linear Power Spectrum
    pyplot.loglog(k_array, h.linear_power(k_array)*k_array**3/(4*numpy.pi),
                  'b--', label='Linear')
    ### Matter-Matter Non-Linear Power Spectrum
    pyplot.loglog(k_array, h.power_mm(k_array)*k_array**3/(4*numpy.pi),
                  'r', label='Matter-Matter')
    ### Galaxy-Matter cross power spectrum using the Zheng HOD
    pyplot.loglog(k_array, h.power_gm(k_array)*k_array**3/(4*numpy.pi),
                  'g-.', label='Galaxy-Matter')
    ### Galaxy-Galaxy auto power using the Zheng HOD
    pyplot.loglog(k_array, h.power_gg(k_array)*k_array**3/(4*numpy.pi),
                  'm:', label='Galaxy-Galaxy')
    
    ### recompute the redshift of the halo model for the redshifts in the loop
    ### below
    for z in [0.5, 1.0, 2.0, 3.0]:
        ### set the halo model redshif to z
        h.set_redshift(z)
        ### recompute and plot the power spectra as before
        pyplot.loglog(k_array, h.linear_power(k_array)*k_array**3/(4*numpy.pi),
                      'b--')
        pyplot.loglog(k_array, h.power_mm(k_array)*k_array**3/(4*numpy.pi),
                      'r')
        pyplot.loglog(k_array, h.power_gm(k_array)*k_array**3/(4*numpy.pi),
                      'g-.')
        pyplot.loglog(k_array, h.power_gg(k_array)*k_array**3/(4*numpy.pi),
                      'm:')
    pyplot.xlim(10**-2, 100)
    pyplot.ylim(10**-2, 10**4)
    pyplot.legend(loc=0)
    pyplot.xlabel('k [h/Mpc]')
    pyplot.ylabel(r'$\Delta^2(k)$')
    pyplot.savefig('test_halo.pdf')

def kernel_unit_test():
    import cosmology
    import kernel
    print "/n*************************"
    print "*                       *"
    print "* Testing kernel Module *"
    print "*                       *"
    print "*************************\n\n"

    print "Testing dNdz"
    print "*************************"

    ### To define a Kernel object we need several things first. We need redshift
    ### distributions as well as the corresponding window functions.
    
    ### initialize a array of redshift and cosmologies for plotting purposes
    z_array = numpy.linspace(0.0, 2.0, 100)
    cosmo = cosmology.MultiEpoch(0.0, 2.0)
    
    ### initilized to galaxy redshift distributions one as a magnitude limited
    ### sample, the other a Guassian with mean z=1.0
    lens_dist = kernel.dNdzMagLim(z_min=0.0, z_max=2.0, a=2, z0=0.3, b=2)
    source_dist = kernel.dNdzGaussian(z_min=0.0, z_max=2.0, z0=1.0, sigma_z=0.2)
    ### normalize the distributions and create PDFs
    lens_dist.normalize()
    source_dist.normalize()

    ### plot
    pyplot.figure(figsize=(10, 8))
    pyplot.plot(z_array, lens_dist.dndz(z_array), label="Lens dN/dz")
    pyplot.plot(z_array, source_dist.dndz(z_array), label="Source dN/dz")
    pyplot.legend(loc=0)
    pyplot.savefig('test_dndz.pdf')

    print "Testing WindowFunction"
    print "*************************"

    ### using the distributions defined above compute the distance weighted
    ### window functions for use in projecting a powerspectrum
    ### Define a galaxy window function
    lens_window = kernel.WindowFunctionGalaxy(redshift_dist=lens_dist)
    ### Define a lensed population of galaxies
    source_window = kernel.WindowFunctionConvergence(redshift_dist=source_dist)
    
    ### plot the window functions as a function of redshift
    chi_array = cosmo.comoving_distance(z_array)
    pyplot.figure(figsize=(10, 8))
    pyplot.plot(z_array, lens_window.window_function(chi_array),
                label="Lens Window")
    pyplot.plot(z_array, source_window.window_function(chi_array),
                label="Source Window(Convergence)")
    pyplot.legend(loc=0)
    pyplot.xlabel('z')
    pyplot.savefig('test_window.pdf')

    print "Testing Kernel"
    print "*************************"
    
    ### Initialize the kernel objects for projecting a power spectrum in z space
    ### Initilize the kernel for galaxy clustering
    k_Auto = kernel.Kernel(ktheta_min=0.001*degToRad*0.001, 
                           ktheta_max=1.0*degToRad*100.0,
                           window_function_a=lens_window, 
                           window_function_b=lens_window)
    ### Kernel computing lensing convergence
    k_Con = kernel.Kernel(0.001*degToRad*0.001, 1.0*degToRad*100.0,
                                  lens_window, source_window)
    ### Kernel for computing galaxy-galaxy lensing
    k_Shear = kernel.GalaxyGalaxyLensingKernel(0.001*degToRad*0.001,
                                               1.0*degToRad*100.0,
                                               lens_window, source_window)

    ### plot results. If the kernels are numericaly noisy for large k*theta
    ### decrease the value (increase the prescision) of the varriable 
    ### kernel_precision in defaults.py
    pyplot.figure(figsize=(10, 8))
    pyplot.semilogx(numpy.exp(k_Auto._ln_ktheta_array)/degToRad,
                    k_Auto.kernel(k_Auto._ln_ktheta_array),
                    label="Autocorrelation Kernel")
    pyplot.semilogx(numpy.exp(k_Con._ln_ktheta_array)/degToRad,
                    k_Con.kernel(k_Con._ln_ktheta_array),
                    label="Convergence Kernel")
    pyplot.semilogx(numpy.exp(k_Shear._ln_ktheta_array)/degToRad,
                    k_Shear.kernel(k_Shear._ln_ktheta_array),
                    label="Shear Kernel")
    pyplot.legend(loc=0)
    pyplot.xlabel(r'$k*\theta$ [h*deg/Mpc]')
    pyplot.savefig('test_kernel.pdf')

    ### Print out the redshifts for which the kernel is maximaly sensitive
    print "Redshifts:", k_Auto.z_bar, k_Con.z_bar,k_Shear.z_bar

def hod_unit_test():
    pass

def correlation_unit_test():
    import correlation
    import hod
    import kernel
    print "/n******************************"
    print "*                            *"
    print "* Testing Correlation Module *"
    print "*                            *"
    print "******************************\n\n"

    ### Definining a correlation object requires first two window functions
    ### (these could in principle be the same window function), theta bounds to
    ### compute the correlation over, and optionaly and HOD object and a 
    ### specification as to which power spectrum to use. Note: the different
    ### correlation classes have approprate default values.

    ### As in kernel_unit_test create galaxy distributions
    lens_dist = kernel.dNdzMagLim(0.0, 2.0, 2, 0.3, 2)
    source_dist = kernel.dNdzGaussian(0.0, 2.0, 1.0, 0.2)

    ### create appropreate window objects
    lens_window = kernel.WindowFunctionGalaxy(lens_dist)
    source_window = kernel.WindowFunctionConvergence(source_dist)

    ### define an hod (optional but needed in order to use power_gm or power_gg)
    zheng = hod.HODZheng(10**13.0, 0.15, 10**13.0, 10**14.0, 1.0)
    ### Define the correlation objects. Note that each of these correlations
    ### is computed using the nonlinear dark matter power spectrum, power_mm.
    ### other options are linear_power which computes the correlation for the
    ### linear spectrum only, power_gm which is the galaxy-matter cross spectrum
    ### and power_gg which is the galaxy-galaxy power spectrum.
    ### Here we define the correlation function of galaxy clustering, note that
    ### it takes only one window function as the second is assumed identical
    auto = correlation.AutoCorrelation(theta_min=0.001*degToRad,
                                       theta_max=1.0*degToRad, 
                                       window_function_galaxy=lens_window,
                                       input_hod=zheng,
                                       powSpec='power_mm')
    ### Define the correlation for galaxy-galaxy magnification. Note it takes
    ### and WindowFunctionGalaxy object and a WindowFunctionConvergence Object
    mag = correlation.MagCorrelation(0.001*degToRad, 1.0*degToRad, 
                                     lens_window, source_window,
                                     input_hod=zheng,
                                     powSpec='power_mm')
    ### Define the correlation for galaxy-galaxy shear. Note it takes
    ### and WindowFunctionGalaxy object and a WindowFunctionConvergence Object
    shear = correlation.GalaxyGalaxyLensing(0.001*degToRad, 1.0*degToRad, 
                                            lens_window, source_window,
                                            input_hod=zheng,
                                            powSpec='power_mm')
    ### Compute the correlation functions between the angular bounds
    auto.compute_correlation()
    mag.compute_correlation()
    shear.compute_correlation()

    ### Plot the dark matter correlations as a function of theta.
    pyplot.figure(figsize=(10, 8))
    pyplot.loglog(auto.theta_array/degToRad, auto.wtheta_array,
                  label="Autocorrelation")
    pyplot.loglog(mag.theta_array/degToRad, mag.wtheta_array/2.0,
                  label="Convergence")
    pyplot.loglog(shear.theta_array/degToRad, shear.wtheta_array,
                  label="Galaxy-Galaxy Lensing")
    pyplot.ylim(10**-4, 10**0)
    pyplot.legend(loc=0)
    pyplot.savefig('test_correlation_matter.pdf')

    ### Reset the power spectra in each correlation object to the approrate
    ### galaxy power spectrum and recompute
    auto.set_power_spectrum('power_gg')
    mag.set_power_spectrum('power_gm')
    shear.set_power_spectrum('power_gm')
    auto.compute_correlation()
    mag.compute_correlation()
    shear.compute_correlation()

    ### Plot the glaxy correlation functions. 
    pyplot.figure(figsize=(10, 8))
    pyplot.loglog(auto.theta_array/degToRad, auto.wtheta_array,
                  label="Autocorrelation")
    pyplot.loglog(mag.theta_array/degToRad, mag.wtheta_array/2.0,
                  label="Convergence")
    pyplot.loglog(shear.theta_array/degToRad, shear.wtheta_array,
                  label="Galaxy-Galaxy Lensing")
    pyplot.ylim(10**-4, 10**0)
    pyplot.legend(loc=0)
    pyplot.savefig('test_correlation_galaxy.pdf')
    
def camb_unit_test():
    pass

if __name__ == "__main__":
    cosmology_unit_test()
    mass_function_unit_test()
    halo_unit_test()
    kernel_unit_test()
    hod_unit_test()
    correlation_unit_test()
    camb_unit_test()
