from Nicola import *
from CaTscript_v3 import *

# V.2 -> Corrected error with Delta RA coordinates

def corrNaN(val, sub=0):  #If value is nan, return sub. Otherwise, return value itself
  if numpy.isnan(val):
    return sub
  else:
    return val


def normSpectra(lambdaSpec, FluxSpec, maskedRegions, D):	#It masks the regions and it fits the spectra with a polynomial of the D degree.
  #Fit spectra continuum
  xFit = []
  yFit = []
  for i in range(len(lambdaSpec)):
    tmpVar = 0
    for j in maskedRegions:
      if (lambdaSpec[i] >= j[0]) and (lambdaSpec[i] <= j[1]):
        tmpVar+=1
    if tmpVar == 0:
      xFit.append(lambdaSpec[i])
      yFit.append(FluxSpec[i])
  coeffPoly=numpy.polyfit(xFit, yFit, D)
  #Build polynomial all over the entire lambdarange
  FittedCurve=numpy.zeros(len(lambdaSpec))
  for i in range(len(lambdaSpec)):
    cont=0
    for j in range(len(coeffPoly)):
      cont+=float(coeffPoly[j])*(numpy.array(lambdaSpec[i])**float(len(coeffPoly)-1-j))
    FittedCurve[i]=cont
  #Dividing spectra per fitted continuum
  FluxSpec_normalized=numpy.array(FluxSpec)/numpy.array(FittedCurve)
  return FluxSpec_normalized, FittedCurve



def checkSpec(manualCheck, spectra, SpecType, normalizeContinuum = False):
  if not(os.path.exists('./'+namegal+'/checkSpec.txt')) and manualCheck:
    if verbose:
      print "Checking spectra"
    plt.ion()
    figure(num=0)
    checkVal=0
    for ii in range(len(spectra.xLambdaRest)):
      if spectra.specFitted[ii] < 0:
        checkVal += 1
    if checkVal > 0:
      clf()	#clean figure
      ax = subplot(111)
      if not(normalizeContinuum):
        if SpecType == 1:
          ax.plot(spectra.xLambdaRest,spectra.specFitted,'b-')
          ax.plot(spectra.xLambdaRest,spectra.yFittedCont,'k--')
          ax.set_ylim([0,1.2*numpy.median(spectra.specFitted)])
        if SpecType == 0:
          ax.plot(spectra.xLambdaRest,spectra.specOriginalClean,'b-')
          ax.plot(spectra.xLambdaRest,spectra.specFitted,'r-')
          ax.plot(spectra.xLambdaRest,spectra.yFittedCont,'k--')
          ax.set_ylim([0,1.2*numpy.median(spectra.specOriginalClean)])
      else:
        if SpecType == 1:
          ax.plot(spectra.xLambdaRest,spectra.specNormalized,'b-')
          ax.plot(spectra.xLambdaRest,spectra.yFittedCont,'k--')
          ax.set_ylim([0,1.2*numpy.median(spectra.specNormalized)])
        if SpecType == 0:
          ax.plot(spectra.xLambdaRest,spectra.specNormalized,'b-')
          ax.plot(spectra.xLambdaRest,spectra.yFittedCont,'k--')
          ax.set_ylim([0,1.2*numpy.median(spectra.specNormalized)])
      ax.set_xlim([8450,8750])
      ax.set_title(spectra.nameObj+' S/N = '+str(spectra.SN))
      answer2 = raw_input("Accept spec? (y/n)")
      if answer2[0] == 'y':
        spectra.check = True
        print "Spectra accepted"
      elif answer2[0] == 'n':
        spectra.check = False
        print "Spectra discarted"
    else:
      stdout.write("\rSpec OK")
      stdout.flush()
      spectra.check = True
    tmpCheck = spectra.check
  return tmpCheck


def relativeCoordinates(RA, Dec, RA_gal, Dec_gal): #In degrees
  # From Huchra+91
  DeltaRA = numpy.sin(numpy.radians(RA-RA_gal))*numpy.cos(numpy.radians(Dec))
  DeltaDec = (numpy.sin(numpy.radians(Dec))*numpy.cos(numpy.radians(Dec_gal))-              numpy.cos(numpy.radians(RA-RA_gal))*numpy.cos(numpy.radians(Dec))*numpy.sin(numpy.radians(Dec_gal)))
  #
  return numpy.degrees(DeltaRA), numpy.degrees(DeltaDec)

###############################################

#Classes

class specJacob:
  def __init__(self, tableSlits, ii, RAGal, DecGal, PA):
    c = 299792.458	#speed light [km/s]
    self.nameObj = tableSlits.field('name')[0]
    #self.RASlit_rel = tableSlits.field('RA')[0]-RAGal
    #self.DecSlit_rel = tableSlits.field('Dec')[0]-DecGal
    self.RASlit_rel, self.DecSlit_rel = relativeCoordinates(tableSlits.field('RA')[0], tableSlits.field('Dec')[0], RAGal, DecGal)
    self.velSlit = tableSlits.field('Vel')[0]
    self.sigmaSlit = tableSlits.field('VelDisp')[0]
    self.errvelSlit = tableSlits.field('errVel')[0]
    self.errsigmaSlit = tableSlits.field('errVelDisp')[0]
    self.specOriginal = tableSlits.field('Flux')[0]
    self.specSky = tableSlits.field('WtdSky')[0]
    self.specFitted = tableSlits.field('BestFit')[0] - tableSlits.field('WtdSky')[0]
    self.specOriginalClean = tableSlits.field('Flux')[0] - tableSlits.field('WtdSky')[0]	#Subtraction sky
    self.xLambda = tableSlits.field('Lambda')[0]
    self.specNormalized = numpy.zeros(len(self.xLambda))
    self.ivar = tableSlits.field('ivar')[0]
    self.var = 1./self.ivar	#Variance
    self.skyivarfudge = tableSlits.field('SkyIvarFudge')[0]
    self.skyindex = tableSlits.field('SkyIndex')[0]
    self.counter = ii
    self.xLambdaRest = self.xLambda*(1.-(self.velSlit/c)) #Doppler effect conversion -> spectra in a rest frame
    self.CaTindex = 0	#Initialized
    self.CaTindexCorr = 0	#Initialized
    self.CaTindexErr = 0	#Initialized
    self.sigmaCorr = 0	#Initialized
    self.distance = 0	#Initialized
    self.SN = tableSlits.field('SN')[0]
    self.xFittedCont = []
    self.yFittedCont = []
    self.filename = tableSlits.field('File')[0]
    self.check = 0	#0 not yet checked, 1 checked and good, 2 checked and bad
    self.MCreal = transpose(tableSlits.field('MC_bestfit')[0])	#It contains the 100 fits of the MC realizations
    self.REFF = numpy.sqrt((self.RASlit_rel)**2+(self.DecSlit_rel)**2)	#DATA (deg)  radial distance from the center
    #Calc elliptical distance
    Xn = self.RASlit_rel*numpy.cos(numpy.pi*(90.-PA0)/180.) - self.DecSlit_rel*numpy.sin(numpy.pi*(90.-PA0)/180.)	#Coordinates respect to the coordinate system defined by the semiaxis of the galaxy
    Yn = self.RASlit_rel*numpy.sin(numpy.pi*(90.-PA0)/180.) + self.DecSlit_rel*numpy.cos(numpy.pi*(90.-PA0)/180.)
    self.RellEFF = numpy.sqrt(b_a*(Xn**2)+((Yn**2)/b_a))  #DATA (deg) radial


class specDF:
  def __init__(self, DF, ii, RAGal, DecGal, PA0, b_a):
    c = 299792.458  #speed light [km/s]
    #import code; code.interact(local=locals())
    try:#if len(shape(DF['NAME'])) == 0:
      self.nameObj = DF['NAME'][ii]
      self.RASlit_rel, self.DecSlit_rel = relativeCoordinates(DF['RA'][ii], DF['DEC'][ii], RAGal, DecGal)
      #self.RASlit_rel, self.DecSlit_rel = DF['RA'][ii]-RAGal, DF['DEC'][ii]-DecGal
      self.velSlit, self.errvelSlit = DF['VEL'][ii], DF['ERRVEL'][ii]
      self.sigmaSlit, self.errsigmaSlit = DF['VELDISP'][ii], DF['ERRVELDISP'][ii]
      #
      self.specOriginal = DF['FLUX'][ii]
      self.specSky = DF['WTDSKY'][ii]
      self.specFitted = DF['BESTFIT'][ii] - DF['WTDSKY'][ii]
      self.specOriginalClean = DF['FLUX'][ii] - DF['WTDSKY'][ii]  #Subtraction sky
      self.xLambda = DF['LAMBDA'][ii]
      self.specNormalized = numpy.zeros(len(self.xLambda))
      self.ivar = DF['IVAR'][ii]
      self.var = 1./self.ivar #Variance
      self.skyivarfudge = DF['SKYIVARFUDGE'][ii]
      self.skyindex = DF['SKYINDEX'][ii]
      self.counter = ii
      self.xLambdaRest = self.xLambda*(1.-(self.velSlit/c)) #Doppler effect conversion -> spectra in a rest frame
      self.CaTindex = 0 #Initialized
      self.CaTindexCorr = 0 #Initialized
      self.CaTindexErr = 0  #Initialized
      self.sigmaCorr = 0  #Initialized
      self.distance = 0 #Initialized
      self.SN = DF['SN'][ii]
      self.xFittedCont = []
      self.yFittedCont = []
      self.filename = DF['FILE'][ii]
      self.check = 0  #0 not yet checked, 1 checked and good, 2 checked and bad
      self.MCreal = transpose(DF['MC_BESTFIT'][ii]) #It contains the 100 fits of the MC realizations
    except:#else:
      self.nameObj = DF['NAME'][ii][0]
      self.RASlit_rel, self.DecSlit_rel = relativeCoordinates(DF['RA'][ii][0], DF['DEC'][ii][0], RAGal, DecGal)
      #self.RASlit_rel, self.DecSlit_rel = DF['RA'][ii][0]-RAGal, DF['DEC'][ii][0]-DecGal
      self.velSlit, self.errvelSlit = DF['VEL'][ii][0], DF['ERRVEL'][ii][0]
      self.sigmaSlit, self.errsigmaSlit = DF['VELDISP'][ii][0], DF['ERRVELDISP'][ii][0]
      #
      self.specOriginal = DF['FLUX'][ii][0]
      self.specSky = DF['WTDSKY'][ii][0]
      self.specFitted = DF['BESTFIT'][ii][0] - DF['WTDSKY'][ii][0]
      self.specOriginalClean = DF['FLUX'][ii][0] - DF['WTDSKY'][ii][0]  #Subtraction sky
      self.xLambda = DF['LAMBDA'][ii][0]
      self.specNormalized = numpy.zeros(len(self.xLambda))
      self.ivar = DF['IVAR'][ii][0]
      self.var = 1./self.ivar #Variance
      self.skyivarfudge = DF['SKYIVARFUDGE'][ii][0]
      self.skyindex = DF['SKYINDEX'][ii][0]
      self.counter = ii
      self.xLambdaRest = self.xLambda*(1.-(self.velSlit/c)) #Doppler effect conversion -> spectra in a rest frame
      self.CaTindex = 0 #Initialized
      self.CaTindexCorr = 0 #Initialized
      self.CaTindexErr = 0  #Initialized
      self.sigmaCorr = 0  #Initialized
      self.distance = 0 #Initialized
      self.SN = DF['SN'][ii][0]
      self.xFittedCont = []
      self.yFittedCont = []
      self.filename = DF['FILE'][ii][0]
      self.check = 0  #0 not yet checked, 1 checked and good, 2 checked and bad
      self.MCreal = transpose(DF['MC_BESTFIT'][ii][0]) #It contains the 100 fits of the MC realizations
    #
    self.REFF = numpy.sqrt((self.RASlit_rel)**2+(self.DecSlit_rel)**2)  #DATA (deg)  radial distance from the center
    #Calc elliptical distance
    Xn = self.RASlit_rel*numpy.cos(numpy.pi*(90.-PA0)/180.) - self.DecSlit_rel*numpy.sin(numpy.pi*(90.-PA0)/180.) #Coordinates respect to the coordinate system defined by the semiaxis of the galaxy
    Yn = self.RASlit_rel*numpy.sin(numpy.pi*(90.-PA0)/180.) + self.DecSlit_rel*numpy.cos(numpy.pi*(90.-PA0)/180.)
    self.RellEFF = numpy.sqrt(b_a*(Xn**2)+((Yn**2)/b_a))  #DATA (deg) radial elliptical distance from the center###############################################



###

def retrieveWavelengthRegions(globalFittingContinuum):
  if verbose: print "A default CaT line selection (Foster 2009) will be used."
  # Rest Wavelength continuum intervals
  WLC1a, WLC1b = 8474., 8483.
  WLC2a, WLC2b =8514., 8526.
  WLC3a, WLC3b = 8563., 8577.
  WLC4a,  WLC4b = 8619., 8642.
  WLC5a, WLC5b = 8680., 8705.
#
  #WLC5a=8695   #Excluding iron line
  #WLC5a=8705   #Excluding entire continuum interval
#
  if globalFittingContinuum: Coont = numpy.array([[WLC1a,WLC1b],[WLC2a,WLC2b],[WLC3a,WLC3b],[WLC4a,WLC4b],[WLC5a,WLC5b]])
  else: Coont = numpy.array([[WLC1a,WLC1b],[WLC2a,WLC2b],[WLC2a,WLC2b],[WLC3a,WLC3b],[WLC4a,WLC4b],[WLC5a,WLC5b]])
#
# Rest Wavelength CaT lines
  WLL1a, WLL1b = 8483., 8513.
  WLL2a, WLL2b = 8527., 8557.
  WLL3a, WLL3b = 8647., 8677.
  Lines = numpy.array([[WLL1a,WLL1b],[WLL2a,WLL2b],[WLL3a,WLL3b]])
  return Coont, Lines


def retrieveSigmaCorrection(Coont, Lines, globalFittingContinuum):
  intrinsicDispersion = 1.5/2.3548  #Intrinsic dispersion of Vazdekis spectra (intrinsic FWHM Vazdekis is 1.5pixels)
  sigmaRange = numpy.array(range(400)[1:])
  if os.path.exists('correctionCaT.dat'):
    answer = raw_input("I found a previous correction file. Do you want to re-find the correction parameters using the template spectra? (y/n)? ")
    if answer[0] == 'y':
      os.system('rm correctionCaT.dat')
      finalCoeff = correctSigma(sigmaRange, Coont, Lines, globalFittingContinuum, intrinsicDispersion)
    else:
      finalCoeff = pickle.load(open('correctionCaT.dat','rb'))
  else:
    finalCoeff = correctSigma(sigmaRange, Coont, Lines, globalFittingContinuum, intrinsicDispersion)
  return finalCoeff

def checkSpectra(listSpec, SpecType, Coont, Lines, globalFittingContinuum, manualCheck):
  #
  #
  coeffT, CheckFile = [], []
# Check pre existing checks
  if os.path.exists('./'+namegal+'/checkSpec.txt'):
    answer = raw_input("Apparently you already checked the spectra. Would you like to do it again (overwriting the previous file)? (y/n)?")
    if answer[0] == 'y':
      os.system('rm ./'+namegal+'/checkSpec.txt')
    else:
      tmpfile = asciidata.open('./'+namegal+'/checkSpec.txt')
      for ii in numpy.arange(len(tmpfile[0])):
        CheckFile.append(bool(tmpfile[0][ii]))
        listSpec[ii].check = bool(tmpfile[0][ii])
#
  Check = []
  for ii in numpy.arange(len(listSpec)):
    if verbose:
      stdout.write("\r Computing spectra  %i / %i " % ((ii)+1, len(listSpec)))
      stdout.flush()
    tmpVAR = True
    #
    while tmpVAR:
#### CaT index calculus
      if SpecType == 1:
        inten = listSpec[ii].specFitted
      else:
        inten = listSpec[ii].specOriginalClean
      #
      if not(globalFittingContinuum):
        coont_fitted = fitCoont2(listSpec[ii].xLambdaRest,inten,Coont,listSpec[ii].ivar)
      elif globalFittingContinuum:
        coont_fitted = fitCoont(listSpec[ii].xLambdaRest,inten,Coont,listSpec[ii].ivar)
      #
      xx = numpy.array(listSpec[ii].xLambdaRest)
      dl = listSpec[ii].xLambdaRest[1]-listSpec[ii].xLambdaRest[0]  #dlambda
      listSpec[ii].CaTindex = CaTindex(listSpec[ii].xLambdaRest,inten,Lines,dl,coont_fitted)
      listSpec[ii].yFittedCont = coont_fitted
      listSpec[ii].xFittedCont = xx
      #
      if os.path.exists('./'+namegal+'/checkSpec.txt'):
        tmpCheck = CheckFile[ii]
        tmpVAR = False
      else:
      #Spec check
        if listSpec[ii].SN < 15:
          tmpCheck = False
        else:
          tmpCheck = checkSpec(manualCheck, listSpec[ii],SpecType,normalizeContinuum=False)
        tmpVAR = False
      Check.append(tmpCheck)
  #
  numpy.savetxt('./'+namegal+'/checkSpec.txt',Check)  #Saving accepted and not accepted spectra
  #
  Check = []
  tmptmp = asciidata.open('./'+namegal+'/checkSpec.txt')
  for ii in tmptmp[0]:
    Check.append(ii)
  print '\n'
  return Check



def savePlotSpectra(listSpec, Coont, Lines, SpecType, EllRAD, EllDecD):
  #
  if (os.path.exists('./'+namegal+'/Spec')):
    os.system('rm -R ./'+namegal+'/Spec')
  os.system('mkdir ./'+namegal+'/Spec')
  #
  check2 = []
  for ii in numpy.arange(len(listSpec)):
    if listSpec[ii].SN > 15: check2.append(ii)
  #
  check2 = numpy.array(check2)
  #
  for ii in range(len(listSpec)):
    if not(os.path.exists('./'+namegal+'/Spec/'+str(ii)+'_'+str(listSpec[ii].CaTindexCorr)[0:6]
         +'_'+str(listSpec[ii].nameObj))):
      dummy = os.system('mkdir ./'+namegal+'/Spec/'+str(ii)+'_'+str(listSpec[ii].CaTindexCorr)[0:6]+'_'+
            str(listSpec[ii].nameObj))
  #
  for ii in range(len(listSpec)):
    if verbose:
      stdout.write("\r\r\rSaving file %i / %i " % (ii+1, len(listSpec)))
      stdout.flush()
#
    plt.ioff()
    fig = figure(figsize=(10,12))
    ax1 = subplot(211)
    ax1.plot(listSpec[ii].xLambdaRest,listSpec[ii].specOriginalClean,'b')
    if SpecType == 1:
      ax1.set_ylim([median(listSpec[ii].specFitted)*0.5,median(listSpec[ii].specFitted)*1.5])
    else:
      ax1.plot(listSpec[ii].xLambdaRest,listSpec[ii].specFitted,'k')
      ax1.set_ylim([median(listSpec[ii].specOriginalClean)*0.5,
        median(listSpec[ii].specOriginalClean)*1.5])
  #Coont fit
    textTitle = (r'$\sigma$ = '+str(listSpec[ii].sigmaSlit)+' SN = '+
                  str(listSpec[ii].SN)+' CaT index = '+str(listSpec[ii].CaTindexCorr))
    ax1.set_xlim([8380,8750])
    suptitle(textTitle)
    grid(True)
    ax2 = subplot(212)
    ax2.plot(listSpec[ii].xLambdaRest,listSpec[ii].specOriginalClean,'k')
    for q in numpy.arange(len(Coont)):
      ax2.axvspan(((Coont)[q][0]),((Coont)[q][1]),facecolor='r',alpha=0.5)
    for q in numpy.arange(len(Lines)):
      ax2.axvspan(((Lines)[q][0]),((Lines)[q][1]),facecolor='g',alpha=0.5)
  #Fitted continuum intervals
    ax2.plot(listSpec[ii].xFittedCont,listSpec[ii].yFittedCont,'k--')
    ax2.set_xlim([8380,8750])
    grid(True)
    if listSpec[ii].check:
      ax2.set_title(listSpec[ii].nameObj+' OK')
    else:
      ax2.set_title(listSpec[ii].nameObj+' Rejected')
    ax2.set_ylabel('Flux')
    ax2.set_xlabel('Wavelenght (Angstrom)')
    savefig('./'+namegal+'/Spec/'+str(ii)+'_'+str(listSpec[ii].CaTindexCorr)[0:6]+'_'+
           str(listSpec[ii].nameObj)+'/spec.eps')
    #
    fileOut = open('./'+namegal+'/Spec/'+str(ii)+'_'+str(listSpec[ii].CaTindexCorr)[0:6]+
                 '_'+str(listSpec[ii].nameObj)+'/slitObj.dat', 'wb')
    pickle.dump(listSpec[ii], fileOut)
    fileOut.close()
  #
    if verbose: print ' DONE!'
#
###SAVING FITS FILE
#
    header = pyfits.Header() #Creation Header file
#Spec
    header.set('simple', 'T')
    header.update('bitpix', -64)
    header.set('naxis',1)
    header.set('ctype', 'linear')
    header.set('crval1', listSpec[ii].xLambdaRest[0], 'first wavelength element')
    header.set('crpix1', 1)
    header.set('cdelt1', listSpec[ii].xLambdaRest[1]-listSpec[ii].xLambdaRest[0], 'step wavelength')
#Obs
    header.set('telescop', 'Keck II')
    header.set('instrume', 'DEIMOS')
    header.set('equinox', 2000.0)
    header.set('object', listSpec[ii].nameObj)
    header.set('file', listSpec[ii].filename)
    header.set('RA_rel', listSpec[ii].RASlit_rel)
    header.set('Dec_rel', listSpec[ii].DecSlit_rel)
    header.set('RA_gal', EllRAD)
    header.set('Dec_gal', EllDecD)
#Kin
    header.set('vel', corrNaN(listSpec[ii].velSlit))
    header.set('errvel', corrNaN(listSpec[ii].errvelSlit))
    header.set('Disp', corrNaN(listSpec[ii].sigmaSlit))
    header.set('errDisp', corrNaN(listSpec[ii].errsigmaSlit))
#SKiMS
    header.set('skyindex', listSpec[ii].skyindex)
    header.set('Cat', listSpec[ii].CaTindex)
    header.set('CatCORR', listSpec[ii].CaTindexCorr)
    header.set('CatERR', listSpec[ii].CaTindexErr)
    header.set('dispcorr', listSpec[ii].sigmaCorr)
    header.set('Relleff', listSpec[ii].RellEFF)
#Other
    header.set('creator', 'Nicola Pastorello')
#
#HEADER SECOND HDU
    headerIVAR = pyfits.Header()
#Spec
    headerIVAR.set('simple', 'T')
    headerIVAR.set('bitpix', -64)
    headerIVAR.set('naxis',1)
    headerIVAR.set('naxis1',4096)
    headerIVAR.set('ctype', 'linear')
    headerIVAR.set('crval1', listSpec[ii].xLambdaRest[0], 'first wavelength element')
    headerIVAR.set('crpix1', 1)
    headerIVAR.set('cdelt1', listSpec[ii].xLambdaRest[1]-listSpec[ii].xLambdaRest[0], 'step wavelength')
#
#Creation HDU data
    primary_hdu = pyfits.PrimaryHDU(data=listSpec[ii].specOriginalClean, header=header)
    ivar_hdu = pyfits.PrimaryHDU(data=listSpec[ii].ivar, header=headerIVAR)
    listH = pyfits.HDUList([primary_hdu], file='spec.fits')
    listH.writeto('./'+namegal+'/Spec/'+str(ii)+'_'+str(listSpec[ii].CaTindexCorr)[0:6]+'_'+
                str(listSpec[ii].nameObj)+'/spec.fits')
    listHivar = pyfits.HDUList([ivar_hdu], file='ivar.fits')
    listHivar.writeto('./'+namegal+'/Spec/'+str(ii)+'_'+str(listSpec[ii].CaTindexCorr)[0:6]+'_'+
                      str(listSpec[ii].nameObj)+'/ivar.fits')
  return True




def calc_CaT_errors(listSpec, globalFittingContinuum, Coont, Lines, sigmaCorrection):
  coeffT = []
  for ii in numpy.arange(len(listSpec)):
    tmpCaTvalues = []
    for jj in numpy.arange(len(transpose(listSpec[ii].MCreal))):
      if verbose:
        stdout.write("\rObj n. %i / %i, Realization n. %i / %i" % (ii+1, len(listSpec),
                    jj+1, len(transpose(listSpec[ii].MCreal))))
        stdout.flush()
      inten = transpose(listSpec[ii].MCreal)[jj]
    # Linear fit of the continuum
      if globalFittingContinuum:
        MCfittedContinuum = fitCoont(listSpec[ii].xLambdaRest,inten,Coont,listSpec[ii].ivar)
      else:
        MCfittedContinuum = fitCoont2(listSpec[ii].xLambdaRest,inten,Coont,listSpec[ii].ivar)
      #
      xx = numpy.array(listSpec[ii].xLambdaRest)
      dl = listSpec[ii].xLambdaRest[1]-listSpec[ii].xLambdaRest[0]  #dlambda
      tmpCaTvalues.append(CaTindex(listSpec[ii].xLambdaRest,inten,Lines,dl,MCfittedContinuum))
    #
    if sigmaCorrection:
      listSpec[ii].CaTindexErr = numpy.std(numpy.array(tmpCaTvalues)*listSpec[ii].sigmaCorr-
                                             listSpec[ii].CaTindexCorr) #Corrected by sigma
    else:
      listSpec[ii].CaTindexErr = numpy.std(numpy.array(tmpCaTvalues)-listSpec[ii].CaTindex) #Not corrected by sigma
    #
  return True



def writeOutputFile(listSpec, sigmaCorrection):
#
  for ii in range(len(listSpec)):
    if (ii == 0): #Initialization files
      fileresults = open('./'+namegal+'/CaTindex.txt','w')
    else:
      fileresults = open('./'+namegal+'/CaTindex.txt','a')
    if sigmaCorrection:
      results = [listSpec[ii].nameObj,listSpec[ii].RASlit_rel,listSpec[ii].DecSlit_rel,
             listSpec[ii].CaTindexCorr, listSpec[ii].CaTindexErr, listSpec[ii].distance,
             listSpec[ii].sigmaSlit,listSpec[ii].errsigmaSlit,listSpec[ii].SN,
             listSpec[ii].skyindex, listSpec[ii].REFF,listSpec[ii].RellEFF, listSpec[ii].check]
    else:
      results = [listSpec[ii].nameObj, listSpec[ii].RASlit_rel,listSpec[ii].DecSlit_rel,
              listSpec[ii].CaTindexCorr, listSpec[ii].CaTindexErr, listSpec[ii].distance,
              listSpec[ii].sigmaSlit,listSpec[ii].errsigmaSlit,listSpec[ii].SN,
              listSpec[ii].skyindex, listSpec[ii].REFF,listSpec[ii].RellEFF, listSpec[ii].check]
  #
    for item in results:
      fileresults.write("%s\t"% item)
    fileresults.write("\n")
    fileresults.close()
#
  #Writing header of the tab
  intes = ('#ID\tRA (deg)\tDec (deg)\tCaTindex\tCaTindexErr\tDistance\tsigma\terrSigma\tSN\tskyIndex\tD Circ (deg)\tD Ellip (deg)\tcheck\n')
  with open('./'+namegal+'/CaTindex.txt', 'r+') as f:
    old = f.read() # read everything in the file
    f.seek(0) # rewind
    f.write(intes + old) # write the new line before
    f.close()
  #
  return True
