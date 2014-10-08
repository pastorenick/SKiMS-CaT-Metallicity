#!/usr/bin/env python
# Filename: radialProfilesAndGradients.py
# Creates the radial metallicity plots and extracts the 
# metallicity gradients within the defined ranges
# Errors are both from the fitting and from MC 


# It reads the OutputMet_corr_CLEAN.txt in the galaxy directory, 
# as well as the extracted profiles' files.
#

# Modified to cut the inner radial profile of NGC1407 at -0.25 dex/dex 
# in order to exclude the metallicity bump from the gradient measure


import glob, pandas, argparse, copy
from Nicola import *
from scipy import stats
from scipy.stats import linregress

#The code remeasures the Dell with the equation
# angleRot = (numpy.pi/180.)*(PA0[ii]-90.)
#  xrot, yrot = (RRA *numpy.cos(angleRot) - DDec * numpy.sin(angleRot), 
#   RRA *numpy.sin(angleRot) + DDec * numpy.cos(angleRot))
# 
# Rell = sqrt( b_a*(xrot**2)+ (yrot**2)/b_a  )
#
def findDell(RA, Dec, PA0, b_a):
  angleRot = (numpy.pi/180.)*(PA0-90.)
  xrot, yrot = (RA *numpy.cos(angleRot) - Dec * numpy.sin(angleRot), 
                RA *numpy.sin(angleRot) + Dec * numpy.cos(angleRot))
  # 
  Rell = numpy.sqrt(b_a*(xrot**2)+(yrot**2)/b_a)
  #
  return Rell


#
def errGradientBS(xx, yy, ww, originalVal, numReal = 1000.):
  # 
  #Shuffling
  #
  import random
  list_xx, list_yy, list_ww = [], [], []
  num_elements = len(xx)
  for ii in numpy.arange(numReal):
    random.seed()
    list_xx_single, list_yy_single, list_ww_single = [], [], []
    for jj in numpy.arange(num_elements):
      indexSel = random.choice(numpy.arange(num_elements))
      list_xx_single.append(xx[indexSel])
      list_yy_single.append(yy[indexSel])
      list_ww_single.append(ww[indexSel])
    list_xx.append(list_xx_single)
    list_yy.append(list_yy_single)
    list_ww.append(list_ww_single)
  #
  #Measuring gradients
  #
  list_gradients = []
  for ii in numpy.arange(int(numReal)):
    if len(list_xx[ii]) > 5:
      coeff, V = numpy.polyfit(list_xx[ii], list_yy[ii], 
                          1, cov=True, w=list_ww[ii])
      list_gradients.append(coeff[0])
   #
    else:
      gradient_inner, intercept, r, prob2, errstd = linregress(list_xx[ii], list_yy[ii])
      list_gradients.append(gradient_inner)
  #
  #Measuring variance
  #
  errm = originalVal-numpy.percentile(list_gradients, 16)
  errp = numpy.percentile(list_gradients, 84) - originalVal 
  return errp, errm




#######################
# v....
# 
#######################

innerRangeLog = [-0.5, 0] #in dex
outerRangeLog = [0, 0.4] #in dex

showPlot = False

#Retrieving galaxy parameters' dictionary
lib_path = os.path.abspath('/Users/npastorello/Desktop/Galaxies/General_Studies')
sys.path.append(lib_path)
from galaxyParametersDictionary_v5 import *


# Main


#Retrieving points and profiles
galnames = glob.glob('NGC*')
dicPathInput_points = {}


for ii in ['NGC4449', 'NGC5907', 'NGC720']:
  galnames.remove(ii)

## POINTS
if len(galnames) == 0:
  print "ERROR, NO INPUT DIRECTORIES FOUND"
else:
  for ii in galnames:
    tmpPath = glob.glob(ii+'/OutputMet_corr_CLEAN.txt')
    dicPathInput_points[ii] = './'+tmpPath[0]
    


## PROFILES
dicPathInput_prof = {}
listMissing = []
if len(galnames) == 0:
  print "ERROR, NO INPUT DIRECTORIES FOUND"
else:
  for ii in galnames:
    try:
      tmpPath = glob.glob(ii+'/Z_radialProfile.txt')
      dicPathInput_prof[ii] = './'+tmpPath[0]
    except:
      listMissing.append(ii)



listGradients = []

# Iteration on single galaxies -> Save data in SKiMS_Z_gradients.dat pickle file
for ii in galnames:
  print ii
# READING DATA POINTS
  fileInput = asciidata.open(dicPathInput_points[ii])
  name, CaT, errCaT = [], [], []
  Dell, SN, RA, Dec = [], [], [], []
  Sigma, errSigma, check = [], [] ,[]
  #
  Z, Z_corr, errpZ, errmZ = [], [], [], []
  #
  for jj in numpy.arange(len(fileInput[0])):
    #
    name.append(fileInput[0][jj])
    CaT.append(float(fileInput[3][jj]))
    errCaT.append(float(fileInput[4][jj]))
    #
#    Dell.append(float(fileInput[7][jj])*3600.)
    SN.append(float(fileInput[8][jj]))
    RA.append(float(fileInput[1][jj])*3600.)
    Dec.append(float(fileInput[2][jj])*3600.)
    Dell.append(findDell(RA[jj], Dec[jj], PA0[ii], b_a[ii]))
    #
    Sigma.append(float(fileInput[5][jj]))
    errSigma.append(float(fileInput[6][jj]))
    #
    Z.append(float(fileInput[9][jj]))
    Z_corr.append(float(fileInput[10][jj]))
    errpZ.append(float(fileInput[11][jj]))
    errmZ.append(float(fileInput[12][jj]))
    if (fileInput[13][jj]) == 'False':
      check.append(False)
    elif (fileInput[13][jj]) == 'True':
      check.append(True)
#
  check = numpy.array(check)
# 
#
# READING RADIAL PROFILE
  fileInput = asciidata.open(dicPathInput_prof[ii])
  R_prof, Z_prof, errZp_prof, errZm_prof = [], [], [], []
  #
  for jj in numpy.arange(len(fileInput[0])):
    #
    R_prof.append(float(fileInput[0][jj]))
    Z_prof.append(float(fileInput[1][jj]))
    errZp_prof.append(float(fileInput[2][jj]))
    errZm_prof.append(float(fileInput[3][jj]))
    #
  R_prof, Z_prof = numpy.array(R_prof), numpy.array(Z_prof)
  errZp_prof, errZm_prof = numpy.array(errZp_prof), numpy.array(errZm_prof)
#
# CLEANING PROFILE FOR NOT RELIABLE CONFIDENCE LIMIT POINTS
#
  selOk = numpy.nonzero((errZp_prof < Z_prof) & (Z_prof < errZm_prof))
  R_prof, Z_prof = R_prof[selOk], Z_prof[selOk]
  errZp_prof, errZm_prof = errZp_prof[selOk], errZm_prof[selOk]
#
#'''
# MODIFY KRIGING RADIAL PROFILES IN ORDER TO HAVE THE SAME SPATIAL 
# EXTENSION OF THE DATAPOINTS
#  Dell = numpy.array(Dell)
#  if ((R_prof[0] < min(Dell[check]) < R_prof[1]) #Kriging 1 - Datapoint 1 - Kriging 2
#      or
#      (min(Dell[check]) < R_prof[0])): #Datapoint 1 - Kriging 1
#    mrel = (Z_prof[1]-Z_prof[0])/(R_prof[1]-R_prof[0])
#    R_prof[0] = min(Dell[check])
#    Z_prof[0] = mrel*R_prof[0]+(Z_prof[1]-mrel*R_prof[1])
#  elif ((R_prof[-2] < max(Dell[check]) < R_prof[-1]) #Kriging -2 - Datapoint -1 - Kriging -1
#      or
#      (max(Dell[check]) > R_prof[-1])):  #Kriging -1 - Datapoint -1
#    mrel = (Z_prof[-2]-Z_prof[-1])/(R_prof[-2]-R_prof[-1])
#    R_prof[-1] = max(Dell[check])
#    Z_prof[-1] = mrel*R_prof[-1]+(Z_prof[-2]-mrel*R_prof[-2])
#'''
# PLOTTING
  if showPlot:
    plt.ion()
  else:
    plt.ioff()
#
  fig = figure(num=0, figsize=(6,10))
  clf()
  #
  title(ii)
  ax = subplot(211)
  ax.set_xlabel(r'R [arcsec]')
  ax.set_ylabel(r'[Z/H] [dex]')
#  '''
#  ax.set_ylabel(r'$\sigma$ [km/s]')
#  '''
  #
  ax.errorbar(numpy.array(Dell)[check], numpy.array(Z_corr)[check], 
                yerr = [numpy.abs(numpy.array(errmZ)[check]),
                        numpy.abs(numpy.array(errpZ)[check])],
                capsize=0, color='k', fmt='o', ms=0)
  ax.scatter(numpy.array(Dell)[check], numpy.array(Z_corr)[check], color='k', label='SKiMS datapoints')
#  '''
#  ax.errorbar(numpy.array(Dell)[check], numpy.array(Sigma)[check], 
#                yerr = numpy.abs(numpy.array(errSigma)[check]),
#                capsize=0, color='k', fmt='o', ms=0)
#  ax.scatter(numpy.array(Dell)[check], numpy.array(Sigma)[check], color='k', label='SKiMS datapoints')
#  '''
#
  ax.plot(R_prof, Z_prof, 'b-', label='Kriging profile')
  ax.fill_between(R_prof, errZm_prof, errZp_prof, color='b', alpha=0.5)
  #
  ax.set_ylim([-1.5, 1.5]), ax.set_xlim([0, numpy.max(numpy.array(Dell)[check])+10])
  ax.plot([Reff[ii], Reff[ii]], [ylim()[0],ylim()[1]], 'k--')
  #
  ax.set_xlabel(r'($R$) [arcsec]')
#
  legend(loc=3)
  #
  ax2 = subplot(212)
  ax2.set_xlabel(r'$\log (R/\rm{R_{e})$ [dex]')
  ax2.set_ylabel(r'[Z/H] [dex]')
  #
  logRad = numpy.log10(numpy.array(Dell)[check]/Reff[ii])
  ax2.errorbar(logRad, numpy.array(Z_corr)[check], 
                yerr = [numpy.abs(numpy.array(errmZ)[check]),
                        numpy.abs(numpy.array(errpZ)[check])],
                capsize=0, color='k', fmt='o', ms=0)
  ax2.scatter(logRad, numpy.array(Z_corr)[check], color='k', label='SKiMS datapoints')
 # '''
 # ax2.errorbar(logRad, numpy.array(Sigma)[check], 
 #               yerr = numpy.array(errSigma)[check],
 #               capsize=0, color='k', fmt='o', ms=0)
 # ax2.scatter(logRad, numpy.array(Sigma)[check], color='k', label='SKiMS datapoints')
 # '''
#
  logRad_prof = numpy.log10(R_prof/Reff[ii])
  ax2.plot(logRad_prof, Z_prof, 'b-', label='Kriging profile')
  ax2.fill_between(logRad_prof, errZm_prof, errZp_prof, color='b', alpha=0.5)
  #
  ax2.set_ylim([-1.5, 1.5]), ax2.set_xlim([-1.2, numpy.max(logRad)+0.1])
  ax2.plot([0, 0], [ylim()[0],ylim()[1]], 'k--')
  #
  ax2.set_xlabel(r'log($R/\rm{R_{e}}$) [dex]')
#
### MEASURING AND STORING GRADIENT ###
  # Fitting the metallicity profile within and without the limits, 
  # just considering the portion of the radial gradients within the two 
  # extreme datapoints (innermost and outermost)
  #
  weights = 1./numpy.sqrt(((errZm_prof-Z_prof)**2.)+((errZp_prof-Z_prof)**2.))
  try:
    #
    minR = numpy.max([numpy.min(logRad), innerRangeLog[0]])
    if ii != 'NGC1407':
      selInner = numpy.nonzero((logRad_prof >= minR) & (logRad_prof < innerRangeLog[1]))
    else: 
      selInner = numpy.nonzero((logRad_prof >= minR) & (logRad_prof < -0.23)) #To exclude the bump
    if len(selInner[0]) > 5:
      coeff, V = numpy.polyfit(logRad_prof[selInner], Z_prof[selInner], 
                          1, cov=True, w=weights[selInner])
      gradient_inner, err_gradient_inner_fit = coeff[0], numpy.sqrt(V[0][0])
   #
    else:
      gradient_inner, intercept, r, prob2, errstd = linregress(logRad_prof[selInner], Z_prof[selInner])
      mx = numpy.mean(logRad_prof[selInner])
      sx2 = numpy.sum((logRad_prof[selInner]-mx)**2)
      err_gradient_inner_fit = errstd * numpy.sqrt(1./sx2)
      #
    err_gradient_inner_BSp, err_gradient_inner_BSm = errGradientBS(logRad_prof[selInner], Z_prof[selInner],
                                         weights[selInner], gradient_inner)
    #
    xx1 = numpy.array([minR, 0])
    ax2.plot(xx1, xx1*coeff[0]+coeff[1], 'r--')
  except:
    "Not possible to measure inner gradient for "+ii
    gradient_inner, err_gradient_inner_fit = nan, nan
    err_gradient_inner_BSp, err_gradient_inner_BSm = nan, nan
    #
  try:
    maxR = numpy.min([numpy.max(logRad), outerRangeLog[1]])
    selOuter = numpy.nonzero((logRad_prof >= outerRangeLog[0]) & (logRad_prof < maxR))
    coeff, V = numpy.polyfit(logRad_prof[selOuter], Z_prof[selOuter], 
                          1, cov=True, w=weights[selOuter])
    gradient_outer, err_gradient_outer_fit = coeff[0], numpy.sqrt(V[0][0])
    err_gradient_outer_BSp, err_gradient_outer_BSm = errGradientBS(logRad_prof[selOuter], 
                                Z_prof[selOuter], weights[selOuter], gradient_outer)
    #
    xx2 = numpy.array([0, maxR])
    ax2.plot(xx2, xx2*coeff[0]+coeff[1], 'r--')
  except:
    "Not possible to measure outer gradient for "+ii
    gradient_outer, err_gradient_outer_fit = nan, nan
    err_gradient_outer_BSp, err_gradient_outer_BSm = nan, nan
    #
    #
  listGradients.append([ii, gradient_inner, err_gradient_inner_fit, 
                                    [err_gradient_inner_BSp, err_gradient_inner_BSm], 
                  gradient_outer, err_gradient_outer_fit,
                                    [err_gradient_outer_BSp, err_gradient_outer_BSm] ])
  #
  fileOut = open('./'+ii+'/SKiMS_Z_gradients.dat', 'wb')
  pickle.dump([ii, gradient_inner, err_gradient_inner_fit, [err_gradient_inner_BSp, err_gradient_inner_BSm],  
                  gradient_outer, err_gradient_outer_fit, [err_gradient_outer_BSp, err_gradient_outer_BSm]], fileOut)
  fileOut.close()
  #
  savefig('./'+ii+'/radial_Z.pdf', bbox_inches='tight')



