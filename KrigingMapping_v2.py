#!/usr/bin/env python
# Filename: KrigingMapping_v2.py
# Run R code to create kriging maps from an input text file and extract the 
# the radial profile.

# It needs KrigingMapping_def.py and the OutputMet_corr_CLEAN.txt in the 
# galaxy directory
#


import glob, copy, argparse
from Nicola import *
from KrigingMapping_def_v2 import *

#######################
# v.1 - it works!
# v.2 - gives the possibility to obtain just CaT, S/N, Z or sigma kriging maps
#######################


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
    theta_CaT, theta_SN, theta_Z, theta_sigma = Theta_Kriging[ii], Theta_Kriging[ii], Theta_Kriging[ii], Theta_Kriging[ii]
##
# Kriging Mapping
##
  if verbose: print "Running Kriging mapping for "+ii
  #
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
         theta_r = theta_Z, coeff_r = 3, savePdf = True, 
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
    #
    dummy = KrigingMapPython('./'+ii+'/Kriging/', ii, genTable_sigma, label='sigma',
                            limits = [0, 250]) #For the visualization
  #
  #
#
# Extracting radial profiles
# 
  if (mode == "CaT") | (mode == "all"):
    linear_prof_RCaT, linear_prof_CaT = radialProfile(ii, './'+ii+'/Kriging/gridKrig_CaT.txt', label='CaT', datapoints = genTable_CaT)
  elif (mode == "SN") | (mode == "all"):
    linear_prof_RSN, linear_prof_SN = radialProfile(ii, './'+ii+'/Kriging/gridKrig_SN.txt', label='SN', datapoints = genTable_SN)
  elif (mode == "Z") | (mode == "all"):
    linear_prof_RZ, linear_prof_Z = radialProfile(ii, './'+ii+'/Kriging/gridKrig_Z.txt', label='Z', datapoints = genTable_Z)
  elif (mode == "sigma") | (mode == "all"):
     linear_prof_Rsigma, linear_prof_sigma = radialProfile(ii, './'+ii+'/Kriging/gridKrig_sigma.txt', label='sigma', datapoints = genTable_sigma)
  #
  # MC errors
  #
  totRealizations = 100
  if (mode == "CaT") | (mode == "all"):
    print "\nFinding MC errors for CaT map"
    radial_errm_CaT, radial_errp_CaT = MCerrors(linear_prof_RCaT, totRealizations, ii, genTable_CaT, 
                   theta_CaT, label='CaT')
  #
  elif (mode == "SN") | (mode == "all"):
    print "\nFinding MC errors for S/N map"
    radial_errm_SN, radial_errp_SN = MCerrors(linear_prof_RSN, totRealizations, ii, genTable_SN, 
                   theta_SN, label='SN')
  #  
  elif (mode == "Z") | (mode == "all"):
    print "\nFinding MC errors for [Z/H] map"
    radial_errm_Z, radial_errp_Z = MCerrors(linear_prof_RZ, totRealizations, ii, genTable_Z, 
                   theta_Z, label='Z')
  #  
  elif (mode == "sigma") | (mode == "all"):
    print "\nFinding MC errors for Sigma map"
    radial_errm_sigma, radial_errp_sigma = MCerrors(linear_prof_Rsigma, totRealizations, ii, genTable_sigma, 
                   theta_sigma, label='sigma')
  #
  # SAVING PROFILES AND ERRORS
  #
  if (mode == "CaT") | (mode == "all"):
    Xradial, Yradial = numpy.array(linear_prof_RCaT), numpy.array(linear_prof_CaT)
    errpYradial, errmYradial = numpy.array(radial_errm_CaT), numpy.array(radial_errp_CaT)
    #
    outTable_CaT = transpose(numpy.array([Xradial, Yradial, errpYradial, errmYradial]))
    #
    fileout = open('./'+ii+'/CaT_radialProfile.txt', 'wb')
    numpy.savetxt(fileout, outTable_CaT, delimiter='\t', header='R (arcsec)\tCaT index (Angstrom)\terrCaT+\terrCaT-')
    fileout.close()
    # 
  #
  elif (mode == "SN") | (mode == "all"):
    Xradial, Yradial = numpy.array(linear_prof_RSN), numpy.array(linear_prof_SN)
    errpYradial, errmYradial = numpy.array(radial_errm_SN), numpy.array(radial_errp_SN)
    #
    outTable_SN = transpose(numpy.array([Xradial, Yradial, errpYradial, errmYradial]))
    #
    fileout = open('./'+ii+'/SN_radialProfile.txt', 'wb')
    numpy.savetxt(fileout, outTable_SN, delimiter='\t', header='R (arcsec)\tS/N\terrSN+\terrSN-')
    fileout.close()
    #
  #
  elif (mode == "Z") | (mode == "all"):
    Xradial, Yradial = numpy.array(linear_prof_RZ), numpy.array(linear_prof_Z)
    errpYradial, errmYradial = numpy.array(radial_errm_Z), numpy.array(radial_errp_Z)
    #
    outTable_Z = transpose(numpy.array([Xradial, Yradial, errpYradial, errmYradial]))
    #
    fileout = open('./'+ii+'/Z_radialProfile.txt', 'wb')
    numpy.savetxt(fileout, outTable_Z, delimiter='\t', header='R (arcsec)\t[Z/H] (dex)\terr[Z/H]+\terr[Z/H]-')
    fileout.close()
    #
  #
  elif (mode == "sigma") | (mode == "all"):
    Xradial, Yradial = numpy.array(linear_prof_Rsigma), numpy.array(linear_prof_sigma)
    errpYradial, errmYradial = numpy.array(radial_errm_sigma), numpy.array(radial_errp_sigma)
    #
    outTable_sigma = transpose(numpy.array([Xradial, Yradial, errpYradial, errmYradial]))
    #
    fileout = open('./'+ii+'/sigma_radialProfile.txt', 'wb')
    numpy.savetxt(fileout, outTable_sigma, delimiter='\t', header='R (arcsec)\tsigma (km/s)\terrsigma+\terrsigma-')
    fileout.close()
    #
  if os.path.exists('./'+ii+'/inProgress'):
    os.remove('./'+ii+'/inProgress')
  # Create flag file (work done)
  open('./'+ii+'/done', 'a').close()
  #
  print "DONE with "+ii+" in "+str(round((time.time() - time0)/60.,2))+" minutes."
