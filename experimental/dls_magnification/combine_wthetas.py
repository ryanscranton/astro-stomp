import numpy
from optparse import OptionParser
from glob import glob
import mle_fitter

def MeanWtheta(data_list):
    sum_wtheta = 0.0
    for data in data_list:
        sum_wtheta += data[:,1]/(1.0*len(data_list))
    return sum_wtheta

def MeanWthetaError(data_list):
    mean_wtheta = MeanWtheta(data_list)
    n_region = 1.0*len(data_list)
    mean_wtheta_error = numpy.zeros(mean_wtheta.shape)
    for data in data_list:
        mean_wtheta_error += (mean_wtheta - data[:,1])**2
    return numpy.sqrt(mean_wtheta_error)*(n_region-1)/(n_region)

def Covariance(data_list):
    n_bins = len(data_list[0])
    n_regions = len(data_list)
    covar_array = numpy.zeros((n_bins,n_bins))
    mean_wtheta = MeanWtheta(data_list)
    for data in data_list:
        for i in xrange(n_bins):
            for j in xrange(n_bins):
                covar_array[i,j] += ((n_regions-1.0)/(n_regions))**2*(
                    (data[i,1]-mean_wtheta[i])*
                    (data[j,1]-mean_wtheta[j]))
    return covar_array

parser = OptionParser()
parser.add_option("--input_tag",dest="input_tag",default="",
                  action="store",type="str",
                  help="Wtheta_tag to load")
parser.add_option("--output_tag",dest="output_tag",default="",
                  action="store",type="str",
                  help="Name appended to output file")
(options, args) = parser.parse_args()

file_name_list = numpy.sort(glob("Wtheta_F?p??_"+options.input_tag)).tolist()
print "# of Files:",len(file_name_list)
data_list = []
for i in xrange(len(file_name_list)):
    gal_gal = numpy.zeros(14)
    gal_rand = numpy.zeros(14)
    rand_gal = numpy.zeros(14)
    rand_rand = numpy.zeros(14)
    hold_name = file_name_list.pop(i)
    for name in file_name_list:
        data = numpy.loadtxt(name)
        gal_gal += data[:,3]
        rand_gal += data[:,4]
        gal_rand += data[:,5]
        rand_rand += data[:,6]
    ang_corr = (gal_gal-rand_gal-gal_rand+rand_rand)/rand_rand
    theta = numpy.loadtxt(file_name_list[0],usecols=(0,))
    data_list.append(numpy.array([theta,ang_corr]).transpose())
    file_name_list.insert(i,hold_name)
theta = data_list[0][:,0]
mean_wtheta = MeanWtheta(data_list)
mean_wtheta_err = MeanWthetaError(data_list)
covar = Covariance(data_list)
print theta
print mean_wtheta
print mean_wtheta_err
#print covar

output = numpy.array([theta,mean_wtheta,mean_wtheta_err]).transpose()
numpy.savetxt("Wtheta_"+options.output_tag,output)
numpy.savetxt("Wcovar_"+options.output_tag,covar)

sig_fit = mle_fitter.MLEOnePar(theta, mean_wtheta, covar)
print "Chi^2 Significance:",sig_fit.chi_squared(0.0)
