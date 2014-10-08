from Nicola import *
from scipy import *

#
# version 3: -> Better python implementation


c = 299792.458	#speed light [km/s]


#Gaussian computed with 10k points and centered in the X interval indentical to the actual spectra
def gaussianCONVOLVE(specY, specX, sigma):
  #gaussian=numpy.exp(-((numpy.array(specX)-median(specX))**2.)/(2.*(sigma**2)))	#Sigma in pixel
  gaussian = (1./(numpy.sqrt(2.*numpy.pi))*sigma)*numpy.exp(-((numpy.array(specX)-median(specX))**2.)/(2.*(sigma**2)))	#Sigma in pixel
  convolvedSpec = numpy.convolve(specY,gaussian,mode='same')
  return convolvedSpec

def sigmaPix(sigmaKM_s,lambdaMEAN,dl, intrinsicDispersion):	#Input: initial sigma, mean lambda in the interval and ratio dpix/dlambda
  #Conversion by doppler effect
  c = 299792.458	#speed light [km/s]
  sigmaPIX = numpy.sqrt(((lambdaMEAN*(sigmaKM_s))/(dl*c))**2.+(intrinsicDispersion)**2)
  return sigmaPIX


def CaTindex(wv,inten,Lines,dl,coontFittedTempl):	#Following Cenarro 2001
  wv = numpy.array(wv)
  CaT = [0,0,0]
  for ii in numpy.arange(len(wv)):
    for jj in numpy.arange(len(Lines)):
      if ((wv[ii] >= Lines[jj][0]) and (wv[ii] < Lines[jj][1])):# and (inten[j] < coontFittedTempl[j])): #The last to avoid a subtraction of the CaT index ->WRONG!
        CaT[jj] += (1.-((inten[ii])/coontFittedTempl[ii]))*dl
  index = 0.4*CaT[0]+CaT[1]+CaT[2]
  return index

def fitCoont(wv,inten,coontrange,ivar):  
#Fit of a polynomial curve with degree=deg for the entire continuum 
#range defined in coontrange. Return the coontinuum flux values in 
#an array with len = len(wv)
  coontFittedArray = numpy.zeros(len(wv))
  #Identifying positional indexes for the enclosed regions:
  index = []
  for ii in numpy.arange(len(wv)):
    for jj in numpy.arange(len(coontrange)):
      if ((wv[ii] >= coontrange[jj][0]) and (wv[ii] < coontrange[jj][1])):
        index.append(ii)
  ###Calc several coefficents:
  eps1, eps2, eps3, eps4, eps5 = 0., 0., 0., 0., 0.
  #
  for ii in index:
    eps1 += ivar[ii]
    eps2 += (wv[ii]*ivar[ii])
    eps3 += ((wv[ii]**2.)*ivar[ii])
    eps4 += (inten[ii]*ivar[ii])
    eps5 += ((inten[ii]*wv[ii])*ivar[ii])
  delta = (eps1*eps3)-(eps2**2.)
  if delta > 0:
    alpha1 = (1./delta)*((eps3*eps4)-(eps2*eps5))
    alpha2 = (1./delta)*((eps1*eps5)-(eps2*eps4))
  else:
    alpha1, alpha2 = 0., 0.
  for ii in numpy.arange(len(wv)):
    coontFittedArray[ii] = alpha1+alpha2*wv[ii]
  return coontFittedArray

def fitCoont2(wv,inten,coontrange,ivar): 
#Fit different polynomial curves with degree=deg 
#for the different coontinuum ranges defined in coontrange. 
#Return the coontinuum flux values in an array with len = len(wv). 
#Where the fitted continuum sections intersecated it uses the mean 
#->GlobalFittingContinuum=0
  coontFittedArray = numpy.zeros(len(wv))
  coontFittedArrayTMP = numpy.zeros((3,len(wv)))
  for ii in numpy.arange(3):	#For the three continuum ranges:
    #Identifying positional indexes for the enclosed regions:
    index = []
    for jj in numpy.arange(len(wv)):
      if (((wv[jj] >= coontrange[ii*2][0]) and (wv[jj] < coontrange[ii*2][1])) 
      	or (((wv[jj] >= coontrange[ii*2+1][0]) and (wv[jj] < coontrange[ii*2+1][1])))):
        index.append(jj)
    ###Calc several coefficents:
    eps1, eps2, eps3, eps4, eps5 = 0., 0., 0., 0., 0.
    for jj in index:
      eps1 += ivar[jj]
      eps2 += (wv[jj]*ivar[jj])
      eps3 += ((wv[jj]**2.)*ivar[jj])
      eps4 += (inten[jj]*ivar[jj])
      eps5 += ((inten[jj]*wv[jj])*ivar[jj])
    delta = (eps1*eps3)-(eps2**2.)
    if delta > 0:
      alpha1 = (1./delta)*((eps3*eps4)-(eps2*eps5))
      alpha2 = (1./delta)*((eps1*eps5)-(eps2*eps4))
    else:
      alpha1, alpha2 = 0., 0.
    for kk in numpy.arange(len(wv)):
      coontFittedArrayTMP[ii][kk]=alpha1+alpha2*wv[kk]	#Contains the three global fitted continuums
  #Where the three continuums overlap I split the interval in two
  for ii in numpy.arange(len(wv)):
    if wv[ii]< mean([coontrange[1][0],coontrange[1][1]]): 
      coontFittedArray[ii] = coontFittedArrayTMP[0][ii]	#Left line
    elif wv[ii]> coontrange[4][0]: 
      coontFittedArray[ii] = coontFittedArrayTMP[2][ii]	#Right line
    else: 
      coontFittedArray[ii] = coontFittedArrayTMP[1][ii]	#Medium line
  return coontFittedArray

def correctSigma(sigmaRange, Coont, Lines, globalFittingContinuum, intrinsicDispersion, pathTempl='./13GyrVazdekisModels/'):
  #Reading Vazdekis spectra
  if verbose: print '\nApplying sigma correction'
  listTempl = asciidata.open(pathTempl+'listTemplate.txt')
  intensTempl, headerTempl, xLambdaTempl, specRes = [], [], [], []
  #
  for ii in range(len(listTempl[0])):
    tmpFile = pyfits.open(pathTempl+listTempl[0][ii])
    intensTempl.append(tmpFile[0].data)
    headerTempl = tmpFile[0].header
    lambdaRange = [(headerTempl['CRVAL1'] - (headerTempl['CRPIX1']-1.)*headerTempl['CDELT1']),(headerTempl['CRVAL1'] - (headerTempl['CRPIX1']-1.)*headerTempl['CDELT1'])+headerTempl['CDELT1']*(headerTempl['NAXIS1']-1.)]
    xLambdaTempl.append(lambdaRange[0]+numpy.array(range(len(tmpFile[0].data)))*headerTempl['CDELT1'])
    specRes.append(headerTempl['CDELT1'])
  #
  # Gaussians creations
  #
  if verbose: print 'Gaussians creations'
  sigmaPx = []
  for ii in sigmaRange:
    sigmaPx.append(sigmaPix(ii,8600.,specRes[0],intrinsicDispersion))
  # Convolution of the spectra
  if verbose: print 'Convolution of the spectra'
  convSpec = []
  for ii in numpy.arange(len(sigmaPx)):
    convSpecTmp = []
    for jj in numpy.arange(len(intensTempl)): 
      tmpg = gaussianCONVOLVE(numpy.array(intensTempl[jj]), xLambdaTempl[jj],
                              sigmaPx[ii]*specRes[0])  #Sigma in Angstrom
      convSpecTmp.append(tmpg)  
    convSpec.append(convSpecTmp)
  #
  # CaT calculus on every template convolved spectra
  #
  if verbose: print 'CaT index calculus'
  #
  if not(os.path.exists('./Plot')): os.system('mkdir Plot')
  #
  meanFunction, coeffT = [], []
  if verbose: print "Fitting polynomial"
  #
  index_sigmaCorr = numpy.zeros([len(intensTempl),len(sigmaPx)])
  for ii in numpy.arange(len(sigmaPx)):
    for jj in numpy.arange(len(intensTempl)):
      if verbose:
      	stdout.write("\r Computing template  %i / %i " % 
      		        ((jj)+(ii*len(intensTempl))+1, len(sigmaPx)*len(intensTempl)))
        stdout.flush()
      inten = convSpec[ii][jj]
  #Linear fit of the continuum
      if globalFittingContinuum: 
      	coontFittedTempl = fitCoont(xLambdaTempl[jj], inten, Coont, 
      		                        numpy.ones(len(xLambdaTempl[jj])))
#      		                       	
      else: 
      	coontFittedTempl = fitCoont2(xLambdaTempl[jj],inten,Coont, 
      		                     numpy.ones(len(xLambdaTempl[jj])))	
      #
      dl = specRes[0] #xLambdaTempl[j][1]-xLambdaTempl[j][0]	#dlambda
      index_sigmaCorr[jj][ii] = CaTindex(xLambdaTempl[jj], inten, Lines, dl, coontFittedTempl)	#Caroline CaT index
#CaT index calc for the original templates
  if verbose: print 'CaT index calculus on original templates'
  #
  index_TemplOrig = numpy.zeros(len(intensTempl))
  for ii in numpy.arange(len(intensTempl)):
    inten = intensTempl[ii]
 #Continuum fit
    if globalFittingContinuum: 
      coontFittedTempl = fitCoont(xLambdaTempl[ii], inten, Coont, 
                          numpy.ones(len(xLambdaTempl[ii])))	#the vazdekis spectra are assumed to have a constant variance
    else: 
      coontFittedTempl = fitCoont2(xLambdaTempl[ii], inten, Coont, 
      	                   numpy.ones(len(xLambdaTempl[ii])))
    #
    dl = xLambdaTempl[ii][1]-xLambdaTempl[ii][0]	#dlambda
    index_TemplOrig[ii] = CaTindex(xLambdaTempl[ii], inten, Lines, dl, coontFittedTempl)	#Caroline CaT index
  #
  if verbose: print 'Creating plot'
  fig = figure()
  axtmp = subplot(111)
  indexRatio, index_TemplOrigTmp = [], []
  for ii in numpy.arange(len(index_TemplOrig)):
    indexRatio.append(index_TemplOrig[ii]/numpy.array(index_sigmaCorr[ii]))
    axtmp.plot(sigmaRange,indexRatio[ii])
    index_TemplOrigTmp.append(indexRatio[ii])
  axtmp.set_xlabel(r'$\sigma$ [km/s]')
  axtmp.set_ylabel(r'CaT/CaT$_\sigma$')
  savefig('./'+namegal+'/Plot/ScatterPlot.pdf')
  #
  if verbose: print 'DONE!'
  meanFunction.append(index_TemplOrigTmp)
  meanTempl, finalCoeff = [], []
  #
  if verbose: print 'Interpolation correction curve'
  for ii in numpy.arange(len(meanFunction)):	#every degree
    totTempl = numpy.zeros(len(meanFunction[ii][0]))
    for jj in numpy.arange(len(meanFunction[ii])):	#every template
      totTempl += meanFunction[ii][jj]    
    meanTempl.append(totTempl/len(meanFunction[ii]))
    # Fitting polynomium 3 degree
    finalCoeff.append(numpy.polyfit(sigmaRange,meanTempl[ii],3))
  #finalCoeff contains the polynomial coefficents
#Saving
  fileOut = open('./correctionCaT.dat','wb')
  pickle.dump(finalCoeff, fileOut)
  fileOut.close()
  return finalCoeff