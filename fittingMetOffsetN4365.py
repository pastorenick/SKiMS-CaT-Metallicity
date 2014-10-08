#!/usr/bin/env python
#Filename: FittingMetOffsetN4365.py 
from Nicola import *
import random, glob
from random import choice #For bootstrapping
from scipy import optimize
from scipy.optimize import minimize

#Retrieve dictionary
lib_path = os.path.abspath('../../General_Studies')
sys.path.append(lib_path)
from galaxyParametersDictionary_v5 import *


#############
# FUNCTIONS #
#############
def retrieveNameGalFromPath(stringPath, sepSymbol=''):
  if 'NGC' in stringPath:
    pos0, flag = 0, True
    while flag:
      if stringPath[pos0:pos0+3] == 'NGC':
        flag = False  
      else: 
        pos0 += 1
      if pos0 > len(stringPath):
        flag = False
    #
    for ii in numpy.arange(len(stringPath[pos0:])):
      if sepSymbol:
        if (stringPath[pos0+ii] == sepSymbol):
          pos1 = pos0+ii 
          flag = True
          return stringPath[pos0:pos1]
      else:
        if ((stringPath[pos0+ii] == '/') 
          or (stringPath[pos0+ii] == '_')):
          pos1 = pos0+ii 
          flag = True
          return stringPath[pos0:pos1]
    #
    if not(flag):
      print 'Galaxy name not found'
      return ''
    #
  else:
    print 'Galaxy name not found'
    return ''

def findDell(RA, Dec, PA0, b_a):
  angleRot = (numpy.pi/180.)*(PA0-90.)
  xrot, yrot = (RA *numpy.cos(angleRot) - Dec * numpy.sin(angleRot), 
                RA *numpy.sin(angleRot) + Dec * numpy.cos(angleRot))
  # 
  Rell = numpy.sqrt(b_a*(xrot**2)+(yrot**2)/b_a)
  #
  return Rell



def minOffset(offset, ySkims, ySauron, errypSkims, errymSkims, errySauron):
  ySkims = numpy.array(ySkims)+offset
  chi_sq = numpy.sum((numpy.abs(ySkims-numpy.array(ySauron)[:,0]))**2./(numpy.sqrt((errypSkims**2.)+(errymSkims**2.)+(numpy.array(errySauron)[:,0]**2.))))
  return chi_sq

def bootstrapError(ySel, ySauronSel, errypSel, errymSel, errySauronSel, offsetreal, nreal = 100.):  #returns bootstrapping error on offset evaluation
  ll = len(ySel)
  listIndices = numpy.arange(ll)
  listOffsets = []
  for ii in range(int(nreal)):
    newlist = []
    random.seed()
    for ii in range(ll):
      newlist.append(choice(listIndices))
    #Measuring offset
    try:
      offset = scipy.optimize.minimize(minOffset, offsetreal, 
      args=(ySel[newlist], numpy.array(ySauronSel)[newlist], errypSel[newlist], errymSel[newlist], numpy.array(errySauronSel)[newlist]), method='BFGS')
      listOffsets.append(offset.x)
    except:
      offset = scipy.optimize.fmin_bfgs(minOffset, offsetreal, 
      args=(ySel[newlist], numpy.array(ySauronSel)[newlist], errypSel[newlist], errymSel[newlist], numpy.array(errySauronSel)[newlist]))
      listOffsets.append(offset[0])
  return numpy.std(listOffsets) #standard deviation sigma

#############
# MAIN
#############

ii= 'NGC4365'
# ATLAS3d
#
inputFile = '../MergingGC_and_SKiMS_Z/SAURONdata/NGC4365_pop_prof_Duncan.txt'
ATLASdata = asciidata.open(inputFile)
SAURON_X, SAURON_Y, SAURON_Yerr = [], [], []
for jj in numpy.arange(len(ATLASdata[0])):
  SAURON_X.append(numpy.log10(ATLASdata[0][jj]/Reff[ii]))
  SAURON_Y.append(ATLASdata[2][jj])
  SAURON_Yerr.append(ATLASdata[7][jj])

SAURON_X, SAURON_Y, SAURON_Yerr = numpy.array(SAURON_X), numpy.array(SAURON_Y), numpy.array(SAURON_Yerr)
#
# SKiMS
#
inputFile = asciidata.open('./'+ii+'/OutputMet_corr_CLEAN.txt')
name, CaT, errCaT = [], [], []
Dell, SN, RA, Dec = [], [], [], []
Sigma, errSigma, check = [], [] ,[]
#
Z, Z_corr, errpZ, errmZ = [], [], [], []
#
for jj in numpy.arange(len(inputFile[0])):
  #
  name.append(inputFile[0][jj])
  CaT.append(float(inputFile[3][jj]))
  errCaT.append(float(inputFile[4][jj]))
  #
  SN.append(float(inputFile[8][jj]))
  RA.append(float(inputFile[1][jj])*3600.)
  Dec.append(float(inputFile[2][jj])*3600.)
  Dell.append(findDell(RA[jj], Dec[jj], PA0[ii], b_a[ii]))
  #
  Sigma.append(float(inputFile[5][jj]))
  errSigma.append(float(inputFile[6][jj]))
  #
  Z.append(float(inputFile[9][jj]))
  Z_corr.append(float(inputFile[10][jj]))
  errpZ.append(float(inputFile[11][jj]))
  errmZ.append(float(inputFile[12][jj]))
  if (inputFile[13][jj]) == 'False':
    check.append(False)
  elif (inputFile[13][jj]) == 'True':
    check.append(True)
#

check = numpy.array(check)
# 
xSKiMS, ySKiMS = numpy.log10(numpy.array(Dell)[check]/Reff[ii]), numpy.array(Z)[check]
eypSKiMS, eymSKiMS = numpy.array(errpZ)[check], numpy.array(errmZ)[check]
#
#Selection points overlapping with SAURONS
#
overlapIndex = numpy.nonzero(xSKiMS < numpy.max(SAURON_X))
xSel, ySel = xSKiMS[overlapIndex], ySKiMS[overlapIndex]
errpySel, errmySel = eypSKiMS[overlapIndex], eymSKiMS[overlapIndex]
#Finding closer SAURON values
xSauronSel, ySauronSel, yerrSauronSel = [], [], []
for kk in range(len(xSel)):
  tmpArr = []
  for jj in range(len(SAURON_X)):
    tmpArr.append(numpy.abs(SAURON_X[jj]-xSel[kk]))
  selTemp = numpy.nonzero(numpy.array(tmpArr) == numpy.min(numpy.array(tmpArr)))
  xSauronSel.append((SAURON_X)[selTemp])
  ySauronSel.append((SAURON_Y)[selTemp])
  yerrSauronSel.append((SAURON_Yerr)[selTemp])

#Finding best match
from scipy.optimize import minimize
result = scipy.optimize.minimize(minOffset, 1., args=(ySel, ySauronSel, errpySel, errmySel, yerrSauronSel), 
      method='BFGS')
errboot = bootstrapError(ySel, ySauronSel, errpySel, errmySel, yerrSauronSel, result.x)
#


fileOut = open('N4365_ATLASoffset.txt', 'w')
fileOut.write('NGC4365'+'\t'+str(result.x[0])+'\t'+str(errboot))
fileOut.close()
print result.x, errboot



