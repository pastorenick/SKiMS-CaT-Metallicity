#!/usr/bin/env python
# Filename: CaTindexMeasurement_v2.py
# Retrieve the CaT index from SKiMS new output (pickle files)
# requires CaTscript_v3.py, CaTindexMeasurement_def_v2.py

# Default call: CaTindexMeasurement NGC### -v
#

from CaTscript_v3 import *
from CaTindexMeasurement_def_v2 import *
import glob, pandas, argparse

#######################
# v.1 - it works!
# v.2 - fixed issue with pandas dataframe
#
#######################

# Global parameters

parser = argparse.ArgumentParser(description='Measure the CaT index of SKiMS data')
parser.add_argument('galaxy', nargs=1, help='galaxy name (e.g. "NGC####")')
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Verbose mode')
parser.add_argument('-I', '--interactive', action='store_true', default=False, help='Interactive mode')
parser.add_argument('-f', '--fittedspec', action='store_true', default=False, help='Use fitted spectra')
parser.add_argument('-c', '--splitContinuum', action='store_true', default=False, help='Use split continuum')
parser.add_argument('-e', '--MCerrors', action='store_true', default=True, help='Measure the MC errors')
parser.add_argument('-s', '--sigmaCorrection', action='store_true', default=True, help='Apply the velocity dispersion correction to the CaT index')
parser.add_argument('-m', '--manualCheck', action='store_true', default=True, help='Manual check of spectra')


args = parser.parse_args()

__builtins__.namegal = args.galaxy[0]
__builtins__.verbose = args.verbose

globalFittingContinuum = not(args.splitContinuum)
InteractiveParMod = args.interactive
calcErr = args.MCerrors
manualCheck = args.manualCheck
sigmaCorrection = args.sigmaCorrection

if args.fittedspec:
  SpecType = 1 	#Fitted
else:
  SpecType = 0  #Original

'''
namegal = "NGC7457"
SpecType = 0	#0 original, 1 fitted
InteractiveParMod = False
#
globalFittingContinuum = True
manualCheck = True		#Manual check of spectra?
calcErr = True
sigmaCorrection = True
'''
# Main
#currentdir = os.getcwd()

if not(os.path.exists('./'+namegal)): os.mkdir('./'+namegal)

if not(os.path.exists('./'+namegal+'/Plots')): os.mkdir('./'+namegal+'/Plots')

if not(os.path.exists('./'+namegal+'/Outputs')): os.mkdir('./'+namegal+'/Outputs')

#Retrieving galaxy parameters' dictionary
lib_path = os.path.abspath('/Users/npastorello/Desktop/Galaxies/General_Studies')
sys.path.append(lib_path)
from galaxyParametersDictionary_v6 import *


# Retrieving coordinates
#
EllRA, EllDec = CentreCoordinates[namegal]	#Coordinates galaxy
vel0, effR, b_a, PA0 = HelioVel[namegal], Reff[namegal], b_a[namegal], PA0[namegal]
#
# Converting galaxy centre coordinates
(dummy, dummy, dummy, dummy, EllRAH), (dummy, dummy, dummy, dummy, EllDecD) = convAngCoord(EllRA), convAngCoord(EllDec)
EllRAD = EllRAH*15. #now in degree

inputFile = open(glob.glob('./Input/*'+namegal+'*_dataframe.dat')[0], 'rb')

try:
  dataframe = pickle.load(inputFile)
except:
  dataframe = pandas.read_pickle(glob.glob('./Input/*'+namegal+'*_dataframe.dat')[0])

inputFile.close()

# Converting dataframe in list of objects
listSpec = []
#for ii in numpy.arange(len(dataframe['NAME'])):
for ii in dataframe.index:
  tmpObj = specDF(dataframe, ii, EllRAD, EllDecD, PA0, b_a)
  listSpec.append(tmpObj)


# Retrieve restframe wavelength ranges to use for CaT index measurement
Coont, Lines = retrieveWavelengthRegions(globalFittingContinuum)

# Check spectra and measure CaT
#if manualCheck:
Check = checkSpectra(listSpec, SpecType, Coont, Lines, globalFittingContinuum, manualCheck)

# Retrieve and apply Sigma correction
if sigmaCorrection:
  finalCoeff = retrieveSigmaCorrection(Coont, Lines, globalFittingContinuum)
  #Velocity dispersion correction
  finalCorr=[]
  for ii in numpy.arange(len(listSpec)):
    Xcorr = listSpec[ii].sigmaSlit
    Corr = 0.
    for jj in numpy.arange(len(finalCoeff[0])):
      Corr += float(finalCoeff[0][jj])*(numpy.array(Xcorr)**float(len(finalCoeff[0])-1-jj))
    finalCorr.append(Corr)
    listSpec[ii].sigmaCorr = Corr
    listSpec[ii].CaTindexCorr = listSpec[ii].CaTindex*listSpec[ii].sigmaCorr


# Saving plots of all spectra with S/N > 15
savePlotSpectra(listSpec, Coont, Lines, SpecType, EllRAD, EllDecD)

# CaT index Errors calculus
calc_CaT_errors(listSpec, globalFittingContinuum, Coont, Lines, sigmaCorrection)



# Write output file
writeOutputFile(listSpec, sigmaCorrection)


if verbose: print '\nEND\n'
