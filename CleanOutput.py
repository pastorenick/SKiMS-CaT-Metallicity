#!/usr/bin/env python
# Filename: CleanOutput.py
# Read the CaTindexMet_corr.txt file and clean it from duplicates, low S/N 
# objects and spatially (1D and 2D) outliers.

# It saves the OutputMet_corr_CLEAN.txt in the galaxy directory
#

import glob, pandas, argparse, copy
from Nicola import *
from excludedSlits_dictionary import *
#######################
# v.1 - it works!
# 
#######################

visualCheck = False
duplicatesCheck = True
meterrorsCheck = True #Errors on metallicity must be < 0.5 dex
SNCheck = True  #Cut in S/N < 35
DistanceCheck = True  #Clean from points at linear distances > 5 Reff
extremeCaTCheck = True  #Clean from extreme values of CaT (outside the range 2 < CaT < 9)
excludedCheck = True  #Clean from slits in the dictionary of the manually selected slits (excludedSlits_dictionary.py)


def visualizeData(RA, Dec, Dell, Z_corr, errpZ, errmZ, SN, index):
  plt.ion()
  fig = figure(num = 0, figsize=(12.6,5))
  clf()
  ax1, ax2 = subplot(121), subplot(122)
  ax1.set_ylim(-2,1.5)
  ax1.errorbar(Dell*3600., Z_corr, yerr=[numpy.abs(errmZ), numpy.abs(errpZ)], c='k', capsize=0, fmt='.')
  ax1.plot(Dell[index]*3600., Z_corr[index], 'r.', ms=20)
  ax1.set_xlabel(r'$R_{e}$ [arcsec]'), ax1.set_ylabel(r'[Z/H] [dex]')
  #
  ax2.scatter(RA*3600., Dec*3600., c='k')
  ax2.scatter(0, 0, c='g', marker='x')
  ax2.plot(RA[index]*3600., Dec[index]*3600., c='r', ms=200)
  ax2.set_xlabel(r'$\Delta$ RA [arcsec]'), ax2.set_ylabel(r'$\Delta$ Dec [arcsec]')
  #
  ax1.set_title('SN = '+str(SN[index]))
  return True




#Retrieving galaxy parameters' dictionary
lib_path = os.path.abspath('/Users/npastorello/Desktop/Galaxies/General_Studies')
sys.path.append(lib_path)
from galaxyParametersDictionary_v5 import *


# Main

galnames = glob.glob('NGC*')
dicPathInput = {}

if len(galnames) == 0:
  print "ERROR, NO INPUT DIRECTORIES FOUND"
else:
  for ii in galnames:
    tmpPath = glob.glob(ii+'/OutputMet_corr.txt')
    dicPathInput[ii] = './'+tmpPath[0]


tol = 1./3600. #Tolerance
# Iteration on single galaxies
for ii in galnames:
  print "\n##########"
  print ii
  print "##########\n"
  fileInput = asciidata.open(dicPathInput[ii])
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
    Dell.append(float(fileInput[7][jj]))
    SN.append(float(fileInput[8][jj]))
    RA.append(float(fileInput[1][jj]))
    Dec.append(float(fileInput[2][jj]))
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
# Clean duplicates
# 
# Find duplicates
  if duplicatesCheck:
    print "Cleaning duplicates"
    listDoubles = []
    #
    for jj in numpy.arange(len(RA)):
      tmpList = []
      for kk in numpy.arange(jj+1, len(RA)):
        if (abs(RA[jj]-RA[kk]) < tol) and (abs(Dec[jj]-Dec[kk]) < tol):
#          print "Duplicated slits: '"+name[jj]+"' '"+name[kk]+"'"
#          print jj, kk
#          print "S/N: "+str(SN[jj])+" "+str(SN[kk])
#          print "check: "+str(check[jj])+" "+str(check[kk])
#          print "CaT: "+str(CaT[jj])+" "+str(CaT[kk])
          if len(tmpList) == 0:
            tmpList.append(jj)
          tmpList.append(kk)
      if len(tmpList) > 0:
        listDoubles.append(tmpList)
  # Clean duplicates keeping the one with check positive or the highest S/N
    SN = numpy.array(SN)
    for jj in listDoubles:
      import copy
      indexBad = copy.copy(jj)
      # SELECTION OF POSITIVE CHECKS
      tmpIndices = numpy.array(jj)[numpy.nonzero((numpy.array(check)[jj] == 'True') | 
              (numpy.array(check)[jj] == '1'))]
      # SELECTION OF BEST S/N
      if len(tmpIndices) > 0:
        indexGood = numpy.array(tmpIndices)[numpy.nonzero(SN[tmpIndices] == numpy.max(SN[tmpIndices]))]
        indexBad.remove(indexGood[0])
        for kk in indexBad:
          check[kk] = 'False'
  #
    print "DONE"
  #
  #
  if meterrorsCheck:
    print "Cleaning high errors on metallicity"
    for jj in numpy.arange(len(RA)):
      if (numpy.abs(errpZ[jj])+numpy.abs(errmZ[jj])) > 0.5:
        check[jj] = 'False'
    print "DONE"
    #
  #
  #
  if SNCheck:
    print "Cleaning low S/N spectra"
    for jj in numpy.arange(len(RA)):
      if SN[jj] < 35.:
        check[jj] = 'False'
    print "DONE"
    #
  #
  #
  if DistanceCheck:
    print "Cleaning too far spectra"
    for jj in numpy.arange(len(RA)):
      if numpy.sqrt((RA[jj]*3600.)**2.+(Dec[jj]*3600.)**2.) > 5.*Reff[ii]:
        check[jj] = 'False'
    print "DONE"
    #
  #
  #
  if extremeCaTCheck:
    print "Cleaning extreme CaT index values"
    for jj in numpy.arange(len(CaT)):
      if not( 2 < CaT[jj] < 9):
        check[jj] = 'False'
    print "DONE"
    #
  #
  #
  if excludedCheck:
    print "Cleaning manually excluded values in the dictionary"
    if ii in dicExcludedSlits.keys():
      for jj in numpy.arange(len(dicExcludedSlits[ii])):
        #Check that there aren't multiple references to the same slit
        numMatches = len(numpy.nonzero(numpy.array(name) == dicExcludedSlits[ii][jj][1]))
        if numMatches == 1:
          for kk in numpy.arange(len(name)):
            if name[kk] == dicExcludedSlits[ii][jj][1]:
              check[kk] = 'False'
        elif numMatches > 1:
          print "Multiple recorrences for the manually excluded slit "+dicExcludedSlits[ii][jj]
          print "of "+ii
          raw_input('Press a key to continue...')
        else:
          print ii+"'s slit "+dicExcludedSlits[ii][jj]+" not found. "
  #
  #
  if visualCheck:
    print "Visual check"
    check = numpy.array(check)
    selSN = numpy.nonzero((check == '1') | (check == 'True'))
    for jj in numpy.arange(len(selSN[0])):
      dummy = visualizeData(numpy.array(RA)[selSN], numpy.array(Dec)[selSN], 
                          numpy.array(Dell)[selSN], numpy.array(Z_corr)[selSN], 
                          numpy.array(errpZ)[selSN], numpy.array(errmZ)[selSN], 
                          numpy.array(SN)[selSN], jj)
      answer3 = raw_input('Accept (y/n)? ')
      if answer3[0] == 'n':
        check[selSN][jj] = 'False'
    print "DONE"
  #
  # SAVING
  for jj in numpy.arange(len(RA)):
    if (jj == 0): #Initialization files
      fileresults = open('./'+ii+'/OutputMet_corr_CLEAN.txt','w')
    else: 
      fileresults = open('./'+ii+'/OutputMet_corr_CLEAN.txt','a')
      #
    results = [name[jj], RA[jj], Dec[jj], CaT[jj], errCaT[jj], Sigma[jj], errSigma[jj], 
               Dell[jj], SN[jj], Z[jj], Z_corr[jj], errpZ[jj], errmZ[jj], check[jj]] 
      #
    for item in results:
      fileresults.write("%s\t"% item)
    fileresults.write("\n")
    fileresults.close()
  #
  intes = ('#ID\tRA (deg)\tDec (deg)\tCaTindex\tCaTindexErr\tsigma\terrSigma\tD Ellip (deg)\tSN\tMet\tMetCorr\terr+_Met\terr-_Met\tcheck\n')
  with open('./'+ii+'/OutputMet_corr_CLEAN.txt', 'r+') as f:
    old = f.read() # read everything in the file
    f.seek(0) # rewind
    f.write(intes + old) # write the new line before
    f.close()



print "END"
