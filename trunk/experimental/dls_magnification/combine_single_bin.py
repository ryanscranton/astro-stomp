import numpy
from optparse import OptionParser
from glob import glob
import cosmology
import mle_fitter
from matplotlib import pyplot

G = 4.30117902*10**-9
h0 = .704
c = 299792.458
degToRad = numpy.pi/180

def ratio_prediction(z,z_0=4.0):
    chi = cosmo.comoving_distance(z)
    chi_0 = cosmo.comoving_distance(z_0)
    g_z = chi*(chi_0-chi)/chi_0
    return (h0*100/c)**2*(g_z*(1+z))**2

def MeanRatio(data_array):
    ratio_mean = []
    for data in data_array:
        ratio_mean.append(data.sum()/(1.0*len(data)))
    return numpy.array(ratio_mean)

def MeanRatioError(data_array):
    mean_array = MeanRatio(data_array)
    regions = 1.0*data_array.shape[1]
    print regions
    var_array = numpy.empty(data_array.shape[0])
    print mean_array
    for idx, data in enumerate(data_array):
        var_array[idx] = numpy.sum((data - mean_array[idx])**2)
    return numpy.sqrt(var_array)*(regions-1.0)/regions

def Covariance(data_array):
    n_bins = data_array.shape[0]
    n_regions = data_array.shape[1]*1.0
    covar_array = numpy.empty((n_bins,n_bins))
    mean_wtheta = MeanRatio(data_array)
    print mean_wtheta
    for i in xrange(n_bins):
        for j in xrange(n_bins):
            covar_array[i,j] = numpy.sum((mean_wtheta[i]-data_array[i,:])*
                                         (mean_wtheta[j]-data_array[j,:]))
    return covar_array*((n_regions-1)/n_regions)**2

parser = OptionParser()
parser.add_option("--input_tag",dest="input_tag",default="",
                  action="store",type="str",
                  help="Wtheta_tag to load")
parser.add_option("--output_tag",dest="output_tag",default="test",
                  action="store",type="str",
                  help="Name appended to output file")
parser.add_option("--alpha_norm",dest="alpha_norm",default=1.162,
                  action="store",type="float",help="Mag Normalization")
(options, args) = parser.parse_args()

z_array = numpy.array([0.474,0.530,0.590,0.664,0.750,0.841,0.955,1.112])
z_bins = [[0.43,0.5],[0.5,0.56],[0.56,0.63],[0.63,0.7],
          [0.7,0.79],[0.79,0.88],[0.88,1.0],[1.0,1.18]]

cosmo = cosmology.MultiEpoch(0,5.01)

ratios = []
magnification = []
auto_corr = []
ratios = []
for bin,z in zip(z_bins,z_array):
    print "Wtheta_F?p??_"+options.input_tag+'_z'+str(bin[0])+str(bin[1])
    print "Wtheta_F?p??_auto_"+options.input_tag+'_z'+str(bin[0])+str(bin[1])
    mag_name_list = numpy.sort(glob("Wtheta_F?p??_"+options.input_tag+
                                    '_z'+str(bin[0])+str(bin[1]))).tolist()
    auto_name_list = numpy.sort(glob("Wtheta_F?p??_auto_"+options.input_tag+
                                     '_z'+str(bin[0])+str(bin[1]))).tolist()
    print "# of Files:",len(mag_name_list)
    print "# of Files:",len(auto_name_list)

    ratio = []
    for i in xrange(len(mag_name_list)):
        mag_gal_gal = 0.0
        mag_gal_rand = 0.0
        mag_rand_gal = 0.0
        mag_rand_rand = 0.0
        auto_gal_gal = 0.0
        auto_gal_rand = 0.0
        auto_rand_gal = 0.0
        auto_rand_rand = 0.0
        hold_mag_name = mag_name_list.pop(i)
        hold_auto_name = auto_name_list.pop(i)
        for mag, auto in zip(mag_name_list,auto_name_list):
            mag_data = numpy.loadtxt(mag)
            mag_gal_gal += mag_data[3]
            mag_rand_gal += mag_data[4]
            mag_gal_rand += mag_data[5]
            mag_rand_rand += mag_data[6]
            
            auto_data = numpy.loadtxt(auto)
            auto_gal_gal += auto_data[3]
            auto_rand_gal += auto_data[4]
            auto_gal_rand += auto_data[5]
            auto_rand_rand += auto_data[6]
        tmp_mag = ((mag_gal_gal-mag_gal_rand-mag_rand_gal+mag_rand_rand)/
                   mag_rand_rand)
        tmp_auto = ((auto_gal_gal-auto_gal_rand-auto_rand_gal+auto_rand_rand)/
                    auto_rand_rand)
        ratio.append(tmp_mag**2/tmp_auto)
        mag_name_list.insert(i,hold_mag_name)
        auto_name_list.insert(i,hold_auto_name)
    print numpy.mean(ratio)
    ratios.append(ratio)
ratios = numpy.array(ratios)
mean_ratio = MeanRatio(ratios)
mean_ratio_err = MeanRatioError(ratios)
covar = Covariance(ratios)
print mean_ratio
print mean_ratio_err
print covar

sig_fit = mle_fitter.MLEOnePar(numpy.ones(z_array[:7].shape),
                               mean_ratio[:7], covar[:7,:7])
sig_ans = sig_fit.MLE()
print "Ratio Chi^2=",sig_fit.chi_squared(sig_ans)

predict = []
for z in z_array:
    predict.append(ratio_prediction(z,z_0=4.0))
predict = numpy.array(predict)
dm_pow = numpy.loadtxt(
    '/home/morrison/workdir/LBG/Theory/dm_only_dm_1000.ascii')
model_fit = mle_fitter.MLEOnePar((predict*dm_pow)[:7], mean_ratio[:7],
                                 covar[:7,:7])
fit = model_fit.MLE()
print "Model Chi^2=", model_fit.chi_squared(fit)

print dm_pow

pyplot.figure(figsize=(8,8))
pyplot.errorbar(z_array, mean_ratio, mean_ratio_err,fmt='bo',lw=3)
pyplot.plot(z_array, predict*dm_pow*fit,'k--',lw=2,
            label=r'Model \chi^2='+str(model_fit.chi_squared(fit)))
pyplot.plot(z_array, sig_ans*numpy.ones(z_array.shape),'g:',lw=2,
            label='Flat \chi^2='+str(sig_fit.chi_squared(sig_ans)))
pyplot.ylim(-0.01,0.01)
pyplot.legend(loc=0)
pyplot.savefig('jack_ratios_'+options.input_tag+'.pdf')

