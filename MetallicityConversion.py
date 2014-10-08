#!/usr/bin/env python
# Filename: MetallicityConversion.py
# Read the CaTindex.txt file and convert the CaT indices into [Z/H]
# requires conversionZ_H_2.dat

# It saves the OutputMet.txt in the galaxy directory
#

import glob, pandas, argparse
from Nicola import *
#######################
# v.1 - it works!
# 
#######################


#Retrieving galaxy parameters' dictionary
lib_path = os.path.abspath('/Users/npastorello/Desktop/Galaxies/General_Studies')
sys.path.append(lib_path)
from galaxyParametersDictionary_v5 import *

# Retrieve conversion Metallicity
coeff, xx, yySU, xVaz, yVaz = pickle.load(open('conversionZ_H_2.dat', 'rb'))

# Main

galnames = glob.glob('NGC*')
dicPathInput = {}


if len(galnames) == 0:
  print "ERROR, NO INPUT DIRECTORIES FOUND"
else:
  for ii in galnames:
    tmpPath = glob.glob(ii+'/CaTindex.txt')
    dicPathInput[ii] = './'+tmpPath[0]


answer = raw_input('Do you want to convert the metallicities with Pastorello+2014 relation (y/n)? ')



# Iteration on single galaxies
for ii in galnames:
  print ii
  #RETRIEVE CORRECTION
  if answer[0] == 'y':
    if ii in MetOffsetSAURON.keys():
      corr = MetOffsetSAURON[ii][0]
    else:
      corrPar = pickle.load(open('./sigmaCorrection.dat','rb'))
      corr = corrPar[0]*numpy.log10(Sigma0[ii])+corrPar[2]
  #
  fileInput = asciidata.open(dicPathInput[ii])
  name, CaT, errCaT = [], [], []
  Dell, SN, RA, Dec = [], [], [], []
  Sigma, errSigma, check = [], [] ,[]
  #
  angrad = PA0[ii]*numpy.pi/180.
  #
  Z_HSU, Met_corr, errpZ_HSU, errmZ_HSU = [], [], [], []
  #
  for jj in numpy.arange(len(fileInput[0])):
    name.append(fileInput[0][jj])
    CaT.append(fileInput[3][jj])
    Sigma.append(fileInput[6][jj])
    errCaT.append(fileInput[4][jj])
    #
    RArot = fileInput[1][jj]*cos(angrad) - fileInput[2][jj]*sin(angrad)
    Decrot = fileInput[1][jj]*sin(angrad) + fileInput[2][jj]*cos(angrad)
    Dell.append(sqrt(((RArot*(sqrt(b_a[ii])))**2.)+((Decrot/sqrt(b_a[ii]))**2.)))
    SN.append(fileInput[8][jj])
    RA.append(fileInput[1][jj])
    Dec.append(fileInput[2][jj])
    errSigma.append(fileInput[7][jj])
    check.append(fileInput[12][jj])
    #
    # Conversion of CaT into metallicity
    CaTtmp, errCaTtmp = fileInput[3][jj], fileInput[4][jj]
    Ztmp = coeff[0]*(CaTtmp**2.) + coeff[1]*CaTtmp + coeff[2]
    Z_HSU.append(Ztmp)
    #
    errpZtmp = (coeff[0]*((CaTtmp+errCaTtmp)**2.) + 
                coeff[1]*(CaTtmp+errCaTtmp) + coeff[2])
    errpZ_HSU.append(errpZtmp-Ztmp)
    #
    errmZtmp = (coeff[0]*((CaTtmp-errCaTtmp)**2.) + 
                coeff[1]*(CaTtmp-errCaTtmp) + coeff[2])
    errmZ_HSU.append(Ztmp-errmZtmp)
    # applying correction
    Met_corr.append(Ztmp+corr)
    #
    # Saving
    #
    if (jj == 0): #Initialization files
      fileresults = open('./'+ii+'/OutputMet.txt','w')
    else: 
      fileresults = open('./'+ii+'/OutputMet.txt','a')
    #
    results = [name[jj], RA[jj], Dec[jj], CaT[jj], errCaT[jj], Sigma[jj], errSigma[jj], 
               Dell[jj], SN[jj], Z_HSU[jj], errpZ_HSU[jj], errmZ_HSU[jj], check[jj]] 
    #
    for item in results:
      fileresults.write("%s\t"% item)
    fileresults.write("\n")
    fileresults.close()
    ###
    if answer[0] == 'y':
      if (jj == 0): #Initialization files
        fileresults = open('./'+ii+'/OutputMet_corr.txt','w')
      else: 
        fileresults = open('./'+ii+'/OutputMet_corr.txt','a')
      #
      results = [name[jj], RA[jj], Dec[jj], CaT[jj], errCaT[jj], Sigma[jj], errSigma[jj], 
               Dell[jj], SN[jj], Z_HSU[jj], Met_corr[jj], errpZ_HSU[jj], errmZ_HSU[jj], check[jj]] 
      #
      for item in results:
        fileresults.write("%s\t"% item)
      fileresults.write("\n")
      fileresults.close()
#
  #Writing header of the tab
  intes = ('#ID\tRA (deg)\tDec (deg)\tCaTindex\tCaTindexErr\tsigma\terrSigma\tD Ellip (deg)\tSN\tMet\terr+_Met\terr-_Met\tcheck\n')
  with open('./'+ii+'/OutputMet.txt', 'r+') as f:
    old = f.read() # read everything in the file
    f.seek(0) # rewind
    f.write(intes + old) # write the new line before
    f.close()
  ###
  if answer[0] == 'y':
    intes = ('#ID\tRA (deg)\tDec (deg)\tCaTindex\tCaTindexErr\tsigma\terrSigma\tD Ellip (deg)\tSN\tMet\tMetCorr\terr+_Met\terr-_Met\tcheck\n')
    with open('./'+ii+'/OutputMet_corr.txt', 'r+') as f:
      old = f.read() # read everything in the file
      f.seek(0) # rewind
      f.write(intes + old) # write the new line before
      f.close()
  #


if verbose: print '\nEND\n'
