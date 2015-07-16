#!/usr/bin/env python
# Filename: KrigingMapping_v2.py

# This code creates a 2d map using kriging from an input of spatially        #
# scattered data points. In the example, the points are some random          #
# samplings of the NGC 1023 velocity dispersion field.                       #
# A more complete description of the procedure can be found in Appendix A of #
# Pastorello+(2014) (http://adsabs.harvard.edu/abs/2014MNRAS.442.1003P).     #
#
# The code also extracts azimuthally averaged radial profiles for the input  #
# variable, following the galaxy isophotal shape.                            #


# It needs KrigingMapping_def_v3.py and the OutputMet_corr_CLEAN.txt in the
# galaxy directory
#


#
# In order to run the code, the following python packages are requested:     #
#
# -asciidata
# -collections
# -matplotlib
# -numpy
# -os
# -pickle
# -pyfits
# -pylab
# -pyraf
# -random
# -rpy2
# -scipy
# -sys
# -time
#
#
# Copyright: Nicola Pastorello (2015) <nicola.pastorello@gmail.com>
# Most recent version: 16/07/15
#
###############################################

import glob, copy, argparse
from Nicola import *
from KrigingMapping_def_v4 import *

##############################
#       VERSION HISTORY:
# v.1 - it works!
# v.2 - gives the possibility to obtain just CaT, S/N, Z or sigma kriging maps
# v.3 - creates both the linear and the logarithmic output profiles
# v.4 - retrieves the errors via both the Bootstrapping and the Monte Carlo
#       approaches. The kriging map sizes can now be defined via the
#       'sizePixelMap' variable (default = 80x80 pixels).
#
#
#       CALL EXAMPLE:
#       python KrigingMapping_v4.py galaxy mode -t -f
#
#       INPUT PARAMETERS:
#
#       - galaxy: is the name of the galaxy for which the code will retrieve
#                 the kriging 2D map and radial profile. It has to correspond
#                 to the name of the directory where the input data points'
#                 file is.
#
#       - mode: is the variable for which the kriging map is produced.
#               Possible values are 'Z', 'CaT', 'sigma', 'SN', 'all'.
#
#       - "-f": forces the creation of the map also in case it has already been
#               produced before.
#
#       - "-t": measures a potential 'theta' parameter analytically instead of
#               reading it from a dictionary. Generally it should be around
#               10 arcsec for galaxies < 30Mpc.
#
#
#       ERROR ESTIMATION:
#       The code provides a slightly overestimated uncertainty on the final
#       radial profiles. It uses both a Bootstrapping and a Monte Carlo
#       simulation approach.
#       In the first case it draw from the data point sample with
#       replacement, in the second it creates maps with the same points
#       but values within their errorbars (homogeneously distributed).
#       Such errors are then summed in quadrature (ie overestimated) between
#       them and with the dispersion of values within the same radial bin
#       when azimuthally averaged.
#
#       In the Monte Carlo case, the values on the single points are extracted
#       from a Gaussian that has a sigma equal to the average error. In a
#       future version I plan to create a split Gaussian, with different left
#       and right dispersions. Obviously, this error component is not computed
#       for the S/N maps (S/N data points don't have an error associated).
#
#       The uncertainty on the final radial profiles contains the sum in
#       quadrature of these two error components plus (still in quadrature) the
#       standard deviation of the points within the same radial bin. This
#       latter is not is not completely independent from the other two errors.
#       Therefore, the total final uncertainty is overestimated.
#
###############################


parser = argparse.ArgumentParser(description='Creates kriging maps')
parser.add_argument('galaxy', nargs=1, help='galaxy name (e.g. "NGC####" or "all")')
parser.add_argument('mode', nargs=1, help='"CaT", "SN", "Z", "sigma", "all"')
parser.add_argument('-t', '--theta', action='store_true', default=False,
                     help='measures the average distance between points as range')
parser.add_argument('-f', '--forcing', action='store_true', default=False,
                     help='forcing the mapping even if it has already done')

args = parser.parse_args()

__builtins__.namegal = args.galaxy[0]
__builtins__.mode = args.mode[0]

__builtins__.forcing = args.forcing


######
__builtins__.thetaFromDic = not(args.theta) #Instead of measuring the average distance between the points,
                    # it takes the theta from the dictionary
__builtins__.pathNick = './'
__builtins__.savePDF = True #In queue it doesn't work


# Main

# Reading input file
galnames = glob.glob('NGC*')
dicPathInput = {}

if len(galnames) == 0:
  print "ERROR, NO INPUT DIRECTORIES FOUND"
else:
  for ii in galnames:
    tmpPath = glob.glob(ii+'/OutputMet_corr_CLEAN.txt')
    dicPathInput[ii] = './'+tmpPath[0]


#Removing galaxies already mapped ('done' flag file in directory)
if namegal == 'all':
  for ii in dicPathInput.keys():
    if not(forcing):
      if (len(glob.glob(ii+'/done')) > 0) | (ii in ['NGC4449', 'NGC5907']):
#    if (ii in ['NGC4449', 'NGC5907']):
        dicPathInput.pop(ii, None)
    else:
      if (ii in ['NGC4449', 'NGC5907']):
#    if (ii in ['NGC4449', 'NGC5907']):
        dicPathInput.pop(ii, None)
#
  listGalaxiesToRun = dicPathInput.keys()
else:
  listGalaxiesToRun = [namegal]


# Creating table for kriging (X, Y, Z, errZ for just the positive check elements)
for ii in listGalaxiesToRun:
  time0 = time.time()
  #Reading old file
  print "\n######"
  print ii
  print "######\n"
  # Create flag file (existing during work in progress)
  open('./'+ii+'/inProgress', 'a').close()
  #
  fileInput = asciidata.open(dicPathInput[ii])
  name, CaT, errCaT = [], [], []
  SN, RA, Dec = [], [], []
  #
  Sigma, errSigma = [], []
  #
  Z, Z_corr, errpZ, errmZ = [], [], [], []
  #
  check = []
  #
  for jj in numpy.arange(len(fileInput[0])):
    #
    name.append(fileInput[0][jj])
    CaT.append(float(fileInput[3][jj]))
    errCaT.append(float(fileInput[4][jj]))
    #
    SN.append(float(fileInput[8][jj]))
    RA.append(float(fileInput[1][jj])*3600.)
    Dec.append(float(fileInput[2][jj])*3600.)
    #
    Sigma.append(float(fileInput[5][jj]))
    errSigma.append(float(fileInput[6][jj]))
    #
    Z.append(float(fileInput[9][jj]))
    Z_corr.append(float(fileInput[10][jj]))
    errpZ.append(float(fileInput[11][jj]))
    errmZ.append(float(fileInput[12][jj]))
    check.append(fileInput[13][jj])
  #
  selCheck = numpy.nonzero((numpy.array(check) == '1') | (numpy.array(check) == '1.0') | (numpy.array(check) == 'True'))
  #
  if not(os.path.exists('./'+ii+'/Kriging')):
    os.mkdir('./'+ii+'/Kriging')
  #
  #
  if (mode == "CaT") | (mode == "all"):
    genTable_CaT = transpose(numpy.array([numpy.array(RA)[selCheck], numpy.array(Dec)[selCheck],
             numpy.array(CaT)[selCheck], numpy.array(errCaT)[selCheck]]))
    #Saving new files
    fileout = open('./'+ii+'/Kriging/listElements_CaT.txt', 'wb')
    numpy.savetxt(fileout, genTable_CaT, delimiter='\t', header='x\ty\tz\terrz')
    fileout.close()
    #
  elif (mode == "SN") | (mode == "all"):
    genTable_SN = transpose(numpy.array([numpy.array(RA)[selCheck], numpy.array(Dec)[selCheck],
             numpy.array(SN)[selCheck]]))
    #
    fileout = open('./'+ii+'/Kriging/listElements_SN.txt', 'wb')
    numpy.savetxt(fileout, genTable_SN, delimiter='\t', header='x\ty\tz\terrz')
    fileout.close()
    #
  elif (mode == "Z") | (mode == "all"):
    errZ = numpy.sqrt(numpy.array(errpZ)**2.+numpy.array(errmZ)**2.)
    genTable_Z = transpose(numpy.array([numpy.array(RA)[selCheck], numpy.array(Dec)[selCheck],
             numpy.array(Z_corr)[selCheck], errZ[selCheck]]))
    #
    fileout = open('./'+ii+'/Kriging/listElements_Z.txt', 'wb')
    numpy.savetxt(fileout, genTable_Z, delimiter='\t', header='x\ty\tz\terrz')
    fileout.close()
    #
  elif (mode == "sigma") | (mode == "all"):
    genTable_sigma = transpose(numpy.array([numpy.array(RA)[selCheck], numpy.array(Dec)[selCheck],
             numpy.array(Sigma)[selCheck], numpy.array(errSigma)[selCheck]]))
    #
    fileout = open('./'+ii+'/Kriging/listElements_sigma.txt', 'wb')
    numpy.savetxt(fileout, genTable_sigma, delimiter='\t', header='x\ty\tz\terrz')
    fileout.close()
  #
  if not(thetaFromDic):
    if verbose: print "Retrieve Theta analytically"
  ##
  # Finding average distance between points (weighted by the errors), to
  # define range in kriging's semivariogram
    if (mode == "CaT") | (mode == "all"):
      rangeKriging_CaT = getAverageDistance(genTable_CaT[:,0], genTable_CaT[:,1],
                                        errz = genTable_CaT[:,3])
      theta_CaT = int(rangeKriging_CaT)
    #
    elif (mode == "SN") | (mode == "all"):
      rangeKriging_SN = getAverageDistance(genTable_SN[:,0], genTable_SN[:,1])
      theta_SN = int(rangeKriging_SN)
    #
    elif (mode == "Z") | (mode == "all"):
      rangeKriging_Z = getAverageDistance(genTable_Z[:,0], genTable_Z[:,1],
                                        errz = genTable_Z[:,3])
      theta_Z = int(rangeKriging_Z)
    elif (mode == "sigma") | (mode == "all"):
      rangeKriging_sigma = getAverageDistance(genTable_sigma[:,0], genTable_sigma[:,1],
                                         errz = genTable_sigma[:,3])
      theta_Z = int(rangeKriging_sigma)
  #
  else:
    if verbose: print "Retrieve Theta from dictionary"
    theta_CaT, theta_SN, theta_Z, theta_sigma = Theta_Kriging[ii], Theta_Kriging[ii], Theta_Kriging[ii], Theta_Kriging[ii]
##
# Kriging Mapping
##
  if verbose: print "Running Kriging mapping for "+ii

  if (mode == "CaT") | (mode == "all"):
    dummy = KrigingR('./'+ii+'/Kriging/listElements_CaT.txt', visualize=False,
         theta_r = theta_CaT, coeff_r = 3, savePdf = True,
         pathOutput = './'+ii+'/Kriging/', label='CaT')
    if verbose: print "\t CaT Kriging map done!"
    # Create Kriging map with Python
    dummy = KrigingMapPython('./'+ii+'/Kriging/', ii, genTable_CaT, label='CaT',
                            limits = [3., +8]) #For the visualization
  #
  elif (mode == "SN") | (mode == "all"):
    dummy = KrigingR('./'+ii+'/Kriging/listElements_SN.txt', visualize=False,
         theta_r =  theta_SN, coeff_r = 3, savePdf = True,
         pathOutput = './'+ii+'/Kriging/', label='SN')
    if verbose: print "\t S/N Kriging map done!"
    #
    dummy = KrigingMapPython('./'+ii+'/Kriging/', ii, genTable_SN, label='SN',
                            limits = [35., 100]) #For the visualization
  #
  elif (mode == "Z") | (mode == "all"):
    dummy = KrigingR('./'+ii+'/Kriging/listElements_Z.txt', visualize=False,
         theta_r = theta_Z, coeff_r = 3, savePdf = False,
         pathOutput = './'+ii+'/Kriging/', label='Z')
    if verbose: print "\t Z Kriging map done!"
    #
    dummy = KrigingMapPython('./'+ii+'/Kriging/', ii, genTable_Z, label='Z',
                         limits = [-3., +2]) #For the visualization
  #
  elif (mode == "sigma") | (mode == "all"):
    dummy = KrigingR('./'+ii+'/Kriging/listElements_sigma.txt', visualize=False,
         theta_r = theta_sigma, coeff_r = 3, savePdf = True,
         pathOutput = './'+ii+'/Kriging/', label='sigma')
    if verbose: print "\t Sigma Kriging map done!"

    dummy = KrigingMapPython('./'+ii+'/Kriging/', ii, genTable_sigma, label='sigma',
                            limits = [0, 250]) #For the visualization


# '''
# Extracting radial profiles
# '''

  if (mode == "CaT") | (mode == "all"):
    linear_prof_RCaT, linear_prof_CaT, std_prof_lin_CaT = radialProfileLin(ii, './'+ii+'/Kriging/gridKrig_CaT.txt', label='CaT', datapoints = genTable_CaT)
    log_prof_RCaT, log_prof_CaT, std_prof_log_CaT = radialProfileLog(ii, './'+ii+'/Kriging/gridKrig_CaT.txt', label='CaT', datapoints = genTable_CaT)
  elif (mode == "SN") | (mode == "all"):
    linear_prof_RSN, linear_prof_SN, std_prof_lin_SN = radialProfileLin(ii, './'+ii+'/Kriging/gridKrig_SN.txt', label='SN', datapoints = genTable_SN)
    log_prof_RSN, log_prof_SN, std_prof_log_SN = radialProfileLog(ii, './'+ii+'/Kriging/gridKrig_SN.txt', label='SN', datapoints = genTable_SN)
  elif (mode == "Z") | (mode == "all"):
    linear_prof_RZ, linear_prof_Z, std_prof_lin_Z = radialProfileLin(ii, './'+ii+'/Kriging/gridKrig_Z.txt', label='Z', datapoints = genTable_Z)
    log_prof_RZ, log_prof_Z, std_prof_log_Z = radialProfileLog(ii, './'+ii+'/Kriging/gridKrig_Z.txt', label='Z', datapoints = genTable_Z)
  elif (mode == "sigma") | (mode == "all"):
     linear_prof_Rsigma, linear_prof_sigma, std_prof_lin_sigma = radialProfileLin(ii, './'+ii+'/Kriging/gridKrig_sigma.txt', label='sigma', datapoints = genTable_sigma)
     log_prof_Rsigma, log_prof_sigma, std_prof_log_sigma = radialProfileLog(ii, './'+ii+'/Kriging/gridKrig_sigma.txt', label='sigma', datapoints = genTable_sigma)

# '''
# Computing the error via bootstrapping
# '''

  totRealizations = 1000
  if (mode == "CaT") | (mode == "all"):
    print "\nFinding BS errors for CaT map"
    (BS_radial_lin_errm_CaT, BS_radial_lin_errp_CaT,
      BS_radial_lin_median_CaT, BS_n_elements_lin_CaT,
      BS_radial_log_errm_CaT, BS_radial_log_errp_CaT,
      BS_radial_log_median_CaT, BS_n_elements_log_CaT) = MCerrors(linear_prof_RCaT, log_prof_RCaT,
              totRealizations, ii, genTable_CaT, theta_CaT, label='CaT', mode='BS')
  #
  elif (mode == "SN") | (mode == "all"):
    print "\nFinding BS errors for S/N map"
    (BS_radial_lin_errm_SN, BS_radial_lin_errp_SN,
      BS_radial_lin_median_SN, BS_n_elements_lin_SN,
      BS_radial_log_errm_SN, BS_radial_log_errp_SN,
      BS_radial_log_median_SN, BS_n_elements_log_SN) = MCerrors(linear_prof_RSN, log_prof_RSN,
                              totRealizations, ii, genTable_SN, theta_SN, label='SN', mode='BS')
  #
  elif (mode == "Z") | (mode == "all"):
    print "\nFinding BS errors for [Z/H] map"
    (BS_radial_lin_errm_Z, BS_radial_lin_errp_Z,
      BS_radial_lin_median_Z, BS_n_elements_lin_Z,
      BS_radial_log_errm_Z, BS_radial_log_errp_Z,
      BS_radial_log_median_Z, BS_n_elements_log_Z) = MCerrors(linear_prof_RZ, log_prof_RZ,
                              totRealizations, ii, genTable_Z, theta_Z, label='Z', mode='BS')
  #
  elif (mode == "sigma") | (mode == "all"):
    print "\nFinding BS errors for Sigma map"
    (BS_radial_lin_errm_sigma, BS_radial_lin_errp_sigma,
      BS_radial_lin_median_sigma, BS_n_elements_lin_sigma,
      BS_radial_log_errm_sigma, BS_radial_log_errp_sigma,
      BS_radial_log_median_sigma, BS_n_elements_log_sigma) = MCerrors(linear_prof_Rsigma, log_prof_Rsigma,
                              totRealizations, ii, genTable_sigma, theta_sigma, label='sigma', mode='BS')

# '''
# Computing the error via Monte Carlo simulation on the actual data points values
# '''

  totRealizations = 1000
  if (mode == "CaT") | (mode == "all"):
    print "\nFinding MC errors for CaT map"
    (MC_radial_lin_errm_CaT, MC_radial_lin_errp_CaT,
      MC_radial_lin_median_CaT, MC_n_elements_lin_CaT,
      MC_radial_log_errm_CaT, MC_radial_log_errp_CaT,
      MC_radial_log_median_CaT, MC_n_elements_log_CaT) = MCerrors(linear_prof_RCaT, log_prof_RCaT,
              totRealizations, ii, genTable_CaT, theta_CaT, label='CaT', mode='MC')
  #
  elif (mode == "SN") | (mode == "all"):
    print "\nFinding MC errors for S/N map"
    (MC_radial_lin_errm_SN, MC_radial_lin_errp_SN,
      MC_radial_lin_median_SN, MC_n_elements_lin_SN,
      MC_radial_log_errm_SN, MC_radial_log_errp_SN,
      MC_radial_log_median_SN, MC_n_elements_log_SN) = numpy.ones(8)*NaN
  #
  elif (mode == "Z") | (mode == "all"):
    print "\nFinding MC errors for [Z/H] map"
    (MC_radial_lin_errm_Z, MC_radial_lin_errp_Z,
      MC_radial_lin_median_Z, MC_n_elements_lin_Z,
      MC_radial_log_errm_Z, MC_radial_log_errp_Z,
      MC_radial_log_median_Z, MC_n_elements_log_Z) = MCerrors(linear_prof_RZ, log_prof_RZ,
                              totRealizations, ii, genTable_Z, theta_Z, label='Z', mode='MC')
  #
  elif (mode == "sigma") | (mode == "all"):
    print "\nFinding MC errors for Sigma map"
    (MC_radial_lin_errm_sigma, MC_radial_lin_errp_sigma,
      MC_radial_lin_median_sigma, MC_n_elements_lin_sigma,
      MC_radial_log_errm_sigma, MC_radial_log_errp_sigma,
      MC_radial_log_median_sigma, MC_n_elements_log_sigma) = MCerrors(linear_prof_Rsigma, log_prof_Rsigma,
                              totRealizations, ii, genTable_sigma, theta_sigma, label='sigma', mode='MC')

  #
  # SAVING PROFILES AND ERRORS
  #
  if (mode == "CaT") | (mode == "all"):

    radial_lin_errm_CaT = numpy.sqrt((BS_radial_lin_median_CaT- BS_radial_lin_errm_CaT)**2. + (MC_radial_lin_median_CaT- MC_radial_lin_errm_CaT)**2. + std_prof_lin_CaT**2.)
    radial_lin_errp_CaT = numpy.sqrt((BS_radial_lin_median_CaT- BS_radial_lin_errp_CaT)**2. + (MC_radial_lin_median_CaT- MC_radial_lin_errp_CaT)**2. + numpy.array(std_prof_lin_CaT)**2.)
    radial_log_errm_CaT = numpy.sqrt((BS_radial_log_median_CaT- BS_radial_log_errm_CaT)**2. + (MC_radial_log_median_CaT- MC_radial_log_errm_CaT)**2. + std_prof_log_CaT**2.)
    radial_log_errp_CaT = numpy.sqrt((BS_radial_log_median_CaT- BS_radial_log_errp_CaT)**2. + (MC_radial_log_median_CaT- MC_radial_log_errp_CaT)**2. + numpy.array(std_prof_log_CaT)**2.)

    # LINEAR
    Xradial_lin, Yradial_lin = numpy.array(linear_prof_RCaT), numpy.array(linear_prof_CaT)
    errpYradial_lin, errmYradial_lin = numpy.array(radial_lin_errm_CaT), numpy.array(radial_lin_errp_CaT)
    n_elements_lin = numpy.array(MC_n_elements_lin_CaT)
    #
    outTable_CaT = transpose(numpy.array([Xradial_lin, Yradial_lin, errpYradial_lin, errmYradial_lin, n_elements_lin]))
    #
    fileout = open('./'+ii+'/CaT_radialProfile_lin.txt', 'wb')
    numpy.savetxt(fileout, outTable_CaT, delimiter='\t', header='R (arcsec)\tCaT index (Angstrom)\terrCaT+\terrCaT-')
    fileout.close()
    #
    # LOGARITHMIC
    Xradial_log, Yradial_log = numpy.array(log_prof_RCaT), numpy.array(log_prof_CaT)
    errpYradial_log, errmYradial_log = numpy.array(radial_log_errm_CaT), numpy.array(radial_log_errp_CaT)
    n_elements_log = numpy.array(MC_n_elements_log_CaT)
    #
    outTable_CaT = transpose(numpy.array([Xradial_log, Yradial_log, errpYradial_log, errmYradial_log, n_elements_log]))
    #
    fileout = open('./'+ii+'/CaT_radialProfile_log.txt', 'wb')
    numpy.savetxt(fileout, outTable_CaT, delimiter='\t', header='R (dex)\tCaT index (Angstrom)\terrCaT+\terrCaT-')
    fileout.close()
    #
  #
  elif (mode == "SN") | (mode == "all"):

    radial_lin_errm_SN = numpy.sqrt((BS_radial_lin_median_SN- BS_radial_lin_errm_SN)**2. + numpy.array(std_prof_lin_SN)**2.)
    radial_lin_errp_SN = numpy.sqrt((BS_radial_lin_median_SN- BS_radial_lin_errp_SN)**2. + numpy.array(std_prof_lin_SN)**2.)
    radial_log_errm_SN = numpy.sqrt((BS_radial_log_median_SN- BS_radial_log_errm_SN)**2. + numpy.array(std_prof_log_SN)**2.)
    radial_log_errp_SN = numpy.sqrt((BS_radial_log_median_SN- BS_radial_log_errp_SN)**2. + numpy.array(std_prof_log_SN)**2.)

    # LINEAR
    Xradial_lin, Yradial_lin = numpy.array(linear_prof_RSN), numpy.array(linear_prof_SN)
    errpYradial_lin, errmYradial_lin = numpy.array(radial_lin_errm_SN), numpy.array(radial_lin_errp_SN)
    n_elements_lin = numpy.array(MC_n_elements_lin_SN)
    #
    outTable_SN = transpose(numpy.array([Xradial_lin, Yradial_lin, errpYradial_lin, errmYradial_lin, n_elements_lin]))
    #
    fileout = open('./'+ii+'/SN_radialProfile_lin.txt', 'wb')
    numpy.savetxt(fileout, outTable_SN, delimiter='\t', header='R (dex)\tSN\terrSN+\terrSN-')
    fileout.close()
    #
    # LOGARITHMIC
    Xradial_log, Yradial_log = numpy.array(log_prof_RSN), numpy.array(log_prof_SN)
    errpYradial_log, errmYradial_log = numpy.array(radial_log_errm_SN), numpy.array(radial_log_errp_SN)
    n_elements_log = numpy.array(MC_n_elements_log_SN)
    #
    outTable_SN = transpose(numpy.array([Xradial_log, Yradial_log, errpYradial_log, errmYradial_log, n_elements_log]))
    #
    fileout = open('./'+ii+'/SN_radialProfile_log.txt', 'wb')
    numpy.savetxt(fileout, outTable_SN, delimiter='\t', header='R (dex)\tSN\terrSN+\terrSN-')
    fileout.close()
    #
  #
  elif (mode == "Z") | (mode == "all"):
    radial_lin_errm_Z = numpy.sqrt((BS_radial_lin_median_Z - BS_radial_lin_errm_Z)**2. + (MC_radial_lin_median_Z- MC_radial_lin_errm_Z)**2. + numpy.array(std_prof_lin_Z)**2.)
    radial_lin_errp_Z = numpy.sqrt((BS_radial_lin_median_Z- BS_radial_lin_errp_Z)**2. + (MC_radial_lin_median_Z- MC_radial_lin_errp_Z)**2. + numpy.array(std_prof_lin_Z)**2.)
    radial_log_errm_Z = numpy.sqrt((BS_radial_log_median_Z - BS_radial_log_errm_Z)**2. + (MC_radial_log_median_Z- MC_radial_log_errm_Z)**2. + numpy.array(std_prof_log_Z)**2.)
    radial_log_errp_Z = numpy.sqrt((BS_radial_log_median_Z- BS_radial_log_errp_Z)**2. + (MC_radial_log_median_Z- MC_radial_log_errp_Z)**2. + numpy.array(std_prof_log_Z)**2.)

    # LINEAR
    Xradial_lin, Yradial_lin = numpy.array(linear_prof_RZ), numpy.array(linear_prof_Z)
    errpYradial_lin, errmYradial_lin = numpy.array(radial_lin_errm_Z), numpy.array(radial_lin_errp_Z)
    n_elements_lin = numpy.array(MC_n_elements_lin_Z)
    #
    outTable_Z = transpose(numpy.array([Xradial_lin, Yradial_lin, errpYradial_lin, errmYradial_lin, n_elements_lin]))
    #
    fileout = open('./'+ii+'/Z_radialProfile_lin.txt', 'wb')
    numpy.savetxt(fileout, outTable_Z, delimiter='\t', header='R (dex)\tZ (dex)\terrZ+\terrZ-')
    fileout.close()
    #
    # LOGARITHMIC
    Xradial_log, Yradial_log = numpy.array(log_prof_RZ), numpy.array(log_prof_Z)
    errpYradial_log, errmYradial_log = numpy.array(radial_log_errm_Z), numpy.array(radial_log_errp_Z)
    n_elements_log = numpy.array(MC_n_elements_log_Z)
    #
    outTable_Z = transpose(numpy.array([Xradial_log, Yradial_log, errpYradial_log, errmYradial_log, n_elements_log]))
    #
    fileout = open('./'+ii+'/Z_radialProfile_log.txt', 'wb')
    numpy.savetxt(fileout, outTable_Z, delimiter='\t', header='R (dex)\tZ (dex)\terrZ+\terrZ-')
    fileout.close()
    #
  #
  elif (mode == "sigma") | (mode == "all"):

    radial_lin_errm_sigma = numpy.sqrt((BS_radial_lin_median_sigma- BS_radial_lin_errm_sigma)**2. + (MC_radial_lin_median_sigma- MC_radial_lin_errm_sigma)**2. + numpy.array(std_prof_lin_sigma)**2.)
    radial_lin_errp_sigma = numpy.sqrt((BS_radial_lin_median_sigma- BS_radial_lin_errp_sigma)**2. + (MC_radial_lin_median_sigma- MC_radial_lin_errp_sigma)**2. + numpy.array(std_prof_lin_sigma)**2.)
    radial_log_errm_sigma = numpy.sqrt((BS_radial_log_median_sigma- BS_radial_log_errm_sigma)**2. + (MC_radial_log_median_sigma- MC_radial_log_errm_sigma)**2. + numpy.array(std_prof_log_sigma)**2.)
    radial_log_errp_sigma = numpy.sqrt((BS_radial_log_median_sigma- BS_radial_log_errp_sigma)**2. + (MC_radial_log_median_sigma- MC_radial_log_errp_sigma)**2. + numpy.array(std_prof_log_sigma)**2.)

    # LINEAR
    Xradial_lin, Yradial_lin = numpy.array(linear_prof_Rsigma), numpy.array(linear_prof_sigma)
    errpYradial_lin, errmYradial_lin = numpy.array(radial_lin_errm_sigma), numpy.array(radial_lin_errp_sigma)
    n_elements_lin = numpy.array(MC_n_elements_lin_sigma)
    #
    outTable_sigma = transpose(numpy.array([Xradial_lin, Yradial_lin, errpYradial_lin, errmYradial_lin, n_elements_lin]))
    #
    fileout = open('./'+ii+'/sigma_radialProfile_lin.txt', 'wb')
    numpy.savetxt(fileout, outTable_sigma, delimiter='\t', header='R (arcsec)\tsigma (km/s)\terrsigma+\terrsigma-')
    fileout.close()
    #
    # LOGARITHMIC
    Xradial_log, Yradial_log = numpy.array(log_prof_Rsigma), numpy.array(log_prof_sigma)
    errpYradial_log, errmYradial_log = numpy.array(radial_log_errm_sigma), numpy.array(radial_log_errp_sigma)
    n_elements_log = numpy.array(MC_n_elements_log_sigma)
    #
    outTable_sigma = transpose(numpy.array([Xradial_log, Yradial_log, errpYradial_log, errmYradial_log, n_elements_log]))
    #
    fileout = open('./'+ii+'/sigma_radialProfile_log.txt', 'wb')
    numpy.savetxt(fileout, outTable_sigma, delimiter='\t', header='R (arcsec)\tsigma (km/s)\terrsigma+\terrsigma-')
    fileout.close()
    #
  if os.path.exists('./'+ii+'/inProgress'):
    os.remove('./'+ii+'/inProgress')
  # Create flag file (work done)
  open('./'+ii+'/done', 'a').close()
  #
  print "DONE with "+ii+" in "+str(round((time.time() - time0)/60.,2))+" minutes."
