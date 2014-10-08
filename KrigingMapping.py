#!/usr/bin/env python
# Filename: KrigingMapping.py
# Run R code to create kriging maps from an input text file and extract the 
# the radial profile.

# It needs KrigingMapping_def.py and the OutputMet_corr_CLEAN.txt in the 
# galaxy directory
#

import glob, copy
from Nicola import *
from KrigingMapping_def import *

#######################
# v.1 - it works!
# 
#######################

__builtins__.thetaFromDic = True #Instead of measuring the average distance between the points, 
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
for ii in dicPathInput.keys():
#  if (len(glob.glob(ii+'/done')) > 0) | (ii in ['NGC4449', 'NGC5907']):
  if (ii in ['NGC4449', 'NGC5907']):
    dicPathInput.pop(ii, None)



# Creating table for kriging (X, Y, Z, errZ for just the positive check elements) 
for ii in dicPathInput.keys():
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
    Z.append(float(fileInput[9][jj]))
    Z_corr.append(float(fileInput[10][jj]))
    errpZ.append(float(fileInput[11][jj]))
    errmZ.append(float(fileInput[12][jj]))
    check.append(fileInput[13][jj])
  #
  selCheck = numpy.nonzero((numpy.array(check) == '1') | (numpy.array(check) == '1.0') | (numpy.array(check) == 'True'))
  genTable_CaT = transpose(numpy.array([numpy.array(RA)[selCheck], numpy.array(Dec)[selCheck], 
             numpy.array(CaT)[selCheck], numpy.array(errCaT)[selCheck]]))
  genTable_SN = transpose(numpy.array([numpy.array(RA)[selCheck], numpy.array(Dec)[selCheck], 
             numpy.array(SN)[selCheck]]))
  #
  errZ = numpy.sqrt(numpy.array(errpZ)**2.+numpy.array(errmZ)**2.)
  genTable_Z = transpose(numpy.array([numpy.array(RA)[selCheck], numpy.array(Dec)[selCheck], 
             numpy.array(Z_corr)[selCheck], errZ[selCheck]]))
  #
  #Saving new files
  if not(os.path.exists('./'+ii+'/Kriging')):
    os.mkdir('./'+ii+'/Kriging')
  #
  fileout = open('./'+ii+'/Kriging/listElements_CaT.txt', 'wb')
  numpy.savetxt(fileout, genTable_CaT, delimiter='\t', header='x\ty\tz\terrz')
  fileout.close()
  #
  fileout = open('./'+ii+'/Kriging/listElements_SN.txt', 'wb')
  numpy.savetxt(fileout, genTable_SN, delimiter='\t', header='x\ty\tz\terrz')
  fileout.close()
  #
  fileout = open('./'+ii+'/Kriging/listElements_Z.txt', 'wb')
  numpy.savetxt(fileout, genTable_Z, delimiter='\t', header='x\ty\tz\terrz')
  fileout.close()
  #
  if not(thetaFromDic):
  ##
  # Finding average distance between points (weighted by the errors), to 
  # define range in kriging's semivariogram
    rangeKriging_CaT = getAverageDistance(genTable_CaT[:,0], genTable_CaT[:,1],
                                        errz = genTable_CaT[:,3])
    theta_CaT = int(rangeKriging_CaT)
    #
    rangeKriging_SN = getAverageDistance(genTable_SN[:,0], genTable_SN[:,1])
    theta_SN = int(rangeKriging_SN)
    #
    rangeKriging_Z = getAverageDistance(genTable_Z[:,0], genTable_Z[:,1],
                                        errz = genTable_Z[:,3])
    theta_Z = int(rangeKriging_Z)
  else:
    theta_CaT, theta_SN, theta_Z = Theta_Kriging[ii], Theta_Kriging[ii], Theta_Kriging[ii]
#
# Kriging Mapping
#
  if verbose: print "Running Kriging mapping for "+ii
  dummy = KrigingR('./'+ii+'/Kriging/listElements_CaT.txt', visualize=False, 
         theta_r = theta_CaT, coeff_r = 3, savePdf = True, 
         pathOutput = './'+ii+'/Kriging/', label='CaT')
  if verbose: print "\t CaT Kriging map done!"
  #
  dummy = KrigingR('./'+ii+'/Kriging/listElements_SN.txt', visualize=False, 
         theta_r =  theta_SN, coeff_r = 3, savePdf = True, 
         pathOutput = './'+ii+'/Kriging/', label='SN')
  if verbose: print "\t S/N Kriging map done!"
  #
  dummy = KrigingR('./'+ii+'/Kriging/listElements_Z.txt', visualize=False, 
         theta_r = theta_Z, coeff_r = 3, savePdf = True, 
         pathOutput = './'+ii+'/Kriging/', label='Z')
  if verbose: print "\t Z Kriging map done!"
  # Create Kriging map with Python
  dummy = KrigingMapPython('./'+ii+'/Kriging/', ii, genTable_CaT, label='CaT',
                            limits = [3., +8]) #For the visualization
  #
  dummy = KrigingMapPython('./'+ii+'/Kriging/', ii, genTable_SN, label='SN',
                            limits = [35., 100]) #For the visualization
  #
  dummy = KrigingMapPython('./'+ii+'/Kriging/', ii, genTable_Z, label='Z',
                            limits = [-3., +2]) #For the visualization
#
# Extracting radial profiles
# 
  linear_prof_RCaT, linear_prof_CaT = radialProfile(ii, './'+ii+'/Kriging/gridKrig_CaT.txt', label='CaT')
  linear_prof_RSN, linear_prof_SN = radialProfile(ii, './'+ii+'/Kriging/gridKrig_SN.txt', label='SN')
  linear_prof_RZ, linear_prof_Z = radialProfile(ii, './'+ii+'/Kriging/gridKrig_Z.txt', label='Z')
  #
  # MC errors
  #
  totRealizations = 100
  print "\nFinding MC errors for CaT map"
  radial_errm_CaT, radial_errp_CaT = MCerrors(linear_prof_RCaT, totRealizations, ii, genTable_CaT, 
                   theta_CaT, label='CaT')
  #
  print "\nFinding MC errors for S/N map"
  radial_errm_SN, radial_errp_SN = MCerrors(linear_prof_RSN, totRealizations, ii, genTable_SN, 
                   theta_SN, label='SN')
  #
  print "\nFinding MC errors for [Z/H] map"
  radial_errm_Z, radial_errp_Z = MCerrors(linear_prof_RZ, totRealizations, ii, genTable_Z, 
                   theta_Z, label='Z')
  #
  # SAVING PROFILES AND ERRORS
  #
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
  Xradial, Yradial = numpy.array(linear_prof_RZ), numpy.array(linear_prof_Z)
  errpYradial, errmYradial = numpy.array(radial_errm_Z), numpy.array(radial_errp_Z)
  #
  outTable_Z = transpose(numpy.array([Xradial, Yradial, errpYradial, errmYradial]))
  #
  fileout = open('./'+ii+'/Z_radialProfile.txt', 'wb')
  numpy.savetxt(fileout, outTable_Z, delimiter='\t', header='R (arcsec)\t[Z/H] (dex)\terr[Z/H]+\terr[Z/H]-')
  fileout.close()
  #
  if os.path.exists('./'+ii+'/inProgress'):
    os.remove('./'+ii+'/inProgress')
  # Create flag file (work done)
  open('./'+ii+'/done', 'a').close()
  #
  print "DONE with "+ii+" in "+str(round((time.time() - time0)/60.,2))+" minutes."
