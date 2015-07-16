import pickle, numpy, asciidata, rpy2, os, sys, time, random, collections
import rpy2.robjects as robjects
from Nicola import *
from random import choice	#For bootstrapping
lib_path = os.path.abspath('/Users/npastorello/Desktop/Galaxies/General_Studies')
sys.path.append(lib_path)
from galaxyParametersDictionary_v6 import *


# V.2 uses the new R fields kriging routines
# V.3 has some new functions for the new error estimation
# V.4 map pixel size is modifiable

'''
Functions
'''
# To get the sorting indices of an array
def permutation_indices(data):
     return sorted(numpy.arange(len(data)), key = data.__getitem__)

#To create bootstrap realizations
def bootstrapRealization(genTable, pathOutput, realization): #Input is table to give Kriging
  import random
  lines = []
  for jj in genTable:
    lines.append(jj)
  #
  #Shuffling
  #
  newList = []
  for jj in numpy.arange(len(lines)):
    random.seed()
    select = choice(lines)
    # To avoid duplicates, if the line already exists, the positions RA and Dec are
    # offset by a random value in the range -0.5<D<0.5 arcsec.
    if select in numpy.array(newList):
      select[0] += random.random()-0.5
      select[1] += random.random()-0.5
    #
    if len(select) == 4:
      newList.append([select[0],select[1],select[2],select[3]])
    else:
      newList.append([select[0],select[1],select[2]])
#
  newList = numpy.array(newList)
# Save in dir
  if not(os.path.exists(pathOutput+'/BS'+str(realization))):
    os.mkdir(pathOutput+'/BS'+str(realization))
# Savetxt file
  listTmp = []
  for jj in newList:
    listTmp.append('\t'.join(map(str, jj))) #Join elements of the same line
  fileTMP = open(pathOutput+'/BS'+str(realization)+'/realization_'+str(int(realization))+'_Points.txt', 'wb')
  fileTMP.write("\n".join(listTmp))
  fileTMP.close()
  return True

def MCextraction(genTable, pathOutput, realization): #Create a list of points with the same spatial positions but values within the error range.
# The Gaussian from which the new values are extracted has a sigma = to the actual error (in a future expansion, it will reflect the actual asymmetry of the input uncertainties).
# In case of S/N, this doesn't work.
  import random
  # Create list with as many rows as the number of input points
  lines = []
  for jj in genTable:
    lines.append(jj)
  #
  #Shuffling
  #
  newList = []
  for jj in numpy.arange(len(lines)):
    random.seed()
    realValue, sigmaErr = lines[jj][2], lines[jj][3]
    newValue = numpy.random.normal(loc=realValue, scale=sigmaErr)
    newList.append([lines[jj][0],lines[jj][1], newValue, sigmaErr])
#
  newList = numpy.array(newList)
# Save in dir
  if not(os.path.exists(pathOutput+'/MC'+str(realization))):
    os.mkdir(pathOutput+'/MC'+str(realization))
# Savetxt file
  listTmp = []
  for jj in newList:
    listTmp.append('\t'.join(map(str, jj))) #Join elements of the same line
  fileTMP = open(pathOutput+'/MC'+str(realization)+'/realization_'+str(int(realization))+'_Points.txt', 'wb')
  fileTMP.write("\n".join(listTmp))
  fileTMP.close()
  return True


#Run kriging interpolation in R
def KrigingR(pathInput, theta_r = 10., coeff_r = 3.,
             visualize = False, savePdf = False, verbose = False,
             pathOutput = 'Outputs/', label='', full=False, sizePixelMap=80):
#
  r = robjects.r
  r.library('fields')
# LOADING R FITTING FUNCTION
  fcov_r = robjects.r('''
        fitSV.cov <- function(x1,x2,theta,marginal=FALSE,C=NA){
   # return marginal variance
     if( marginal) { return(rep( 1, nrow( x1)))}

    # find cross covariance matrix
      temp<- exp(-(rdist(x1,x2)/theta)**2)
      if( is.na(C[1])){
          return( temp)}
      else{
          return( temp%*%C)}
      } ''')
#
  robjects.globalenv['coeff'] = coeff_r
  robjects.globalenv['range'] = theta_r
#
  filename = pathInput
  robjects.globalenv['tmp'] = filename
  robjects.globalenv['pathOutput'] = pathOutput
  robjects.globalenv['label'] = label
#
  filetab_r = robjects.r('''filetab <- read.table(tmp)''')
  # Extraction variables from table
  selectiontab_r = robjects.r('''selection <- (filetab$V1 !=0)''')
  x_r = robjects.r('''x <- filetab[selection, "V1"]''')
  y_r = robjects.r('''y <- filetab[selection, "V2"]''')
  z_r = robjects.r('''z <- filetab[selection, "V3"]''')
  #
  if label != 'SN':
    zerr_r = robjects.r('''zerr <- filetab[selection, "V4"]''')
    ww_r = robjects.r('''ww <- ((1./(zerr^2.))/max(1./(zerr^2.)))''')
# Creation Position Matrix
  X_r = robjects.r('''X <- data.matrix(data.frame(x,y))''')
#
# Kriging fit
#  fit_r = robjects.r('''fit <- Krig(X, z, weights = ww, cov.function="fitSV.cov", m=coeff,theta=range)''')
  if label != 'SN':
    fit_r = robjects.r('''fit <- Krig(X, z, cov.function="fitSV.cov",
                        m=coeff, theta=range, weights = ww)''')
  else:
    fit_r = robjects.r('''fit <- Krig(X, z, cov.function="fitSV.cov",
                        m=coeff, theta=range)''')
  if visualize:
    robjects.r(''' summary(fit) ''')
    robjects.r(''' set.panel(2,2) ''')
    robjects.r(''' plot(fit) ''')
#
# Visualization map
  if visualize:
    robjects.r('''set.panel()''')
    zrange_r = robjects.r('''zrange <- c(min(z), max(z))''')
    xrange_r = robjects.r('''xrange <- c(max(x),min(x))''')
    robjects.r('''surface(fit, type="C", xlab='RA [arcsec]', ylab='Dec
                [arcsec]',levels=c(-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1,-1.1,-1.2),
                zlab='[Z/H] [dex]', extrap=FALSE, zlim=zrange, xlim=xrange)''')
    robjects.r('''par(new=T)''')
#POINTS
#color
    numcolors = 128.
    robjects.r('''n.color=64''')
    robjects.r('''pos.data <- z - zrange[1]''')
#
    posdata = robjects.r('''pos.data <- pos.data / (zrange[2]-(zrange[1])) *
                         (n.color-1) + 1''')
    timcol = robjects.r('''tim.colors(n.color)''')
    robjects.r('''cols <- tim.colors(n.color)[pos.data]''')
#
  ##size
    if label != 'SN':
      robjects.r('''sizeRange <- c(min(ww), max(ww))''')
  #I keep the size range within 1 and 5
      robjects.r('''sizes <- ((ww) * 4./(sizeRange[2]-sizeRange[1]))+1''')
      robjects.r('''sizes <- sizes - (max(sizes)-4)+1''')
#
      robjects.r('''points(X, pch=24, bg=cols, cex=sizes)''')
    else:
      robjects.r('''points(X, pch=24, bg=cols)''')
#
  if savePdf:
    robjects.r(''' filename <- paste(pathOutput,"Kriging_",label,".pdf", sep="") ''')
    robjects.r(''' pdf(filename) ''')
    robjects.r('''set.panel()''')
    zrange_r = robjects.r('''zrange <- c(min(z), max(z))''')
    zrange_r = robjects.r('''zrange <- c(-2, 2)''')
    xrange_r = robjects.r('''xrange <- c(max(x),min(x))''')
    robjects.r('''surface(fit, type="C", xlab='RA [arcsec]', ylab='Dec [arcsec]',
                  levels=c(-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1,-1.1,-1.2), zlab='[Z/H]
                  [dex]', extrap=FALSE, zlim=zrange, xlim=xrange)''')
    robjects.r('''par(new=T)''')
    #POINTS
    ##color
    robjects.r('''n.color=64''')
    robjects.r('''pos.data <- z - zrange[1]''')
    robjects.r('''pos.data <- pos.data / (zrange[2]-(zrange[1])) * (n.color-1) + 1''')
    robjects.r('''cols <- tim.colors(n.color)[pos.data]''')
  ##size
    if label != 'SN':
      robjects.r('''sizeRange <- c(min(ww), max(ww))''')
  #I keep the size range within 1 and 5
      robjects.r('''sizes <- ((ww) * 4./(sizeRange[2]-sizeRange[1]))+1''')
      robjects.r('''sizes <- sizes - (max(sizes)-4)+1''')
#
      robjects.r('''points(X, pch=24, bg=cols, cex=sizes)''')
    else:
      robjects.r('''points(X, pch=24, bg=cols)''')
    robjects.r(''' dev.off() ''')
#
# Extration map grid
#
  import platform
  posOutput_surface = 3 #Works only with the old version of 'predict'
  if platform.system() == 'Linux':  #On G2 there is a newer version of fields
    look_r = robjects.r(''' look<-predictSurface(fit,
          nx = '''+str(int(sizePixelMap))+''', ny = '''+str(int(sizePixelMap))+''') ''')
#standardErrorsGrid <- c(predict.se(fit,gridK)) #On a grid
    standardErrorsGrid_r = robjects.r(''' standardErrorsGrid <-
                                    c(predictSurfaceSE(fit,
          nx = '''+str(int(sizePixelMap))+''', ny = '''+str(int(sizePixelMap))+''')) ''')
  else:
    try:
      if full:
        look_r = robjects.r(''' look<-predict.surface(fit, extrap=TRUE,
          nx = '''+str(int(sizePixelMap))+''', ny = '''+str(int(sizePixelMap))+''') ''')
#standardErrorsGrid <- c(predict.se(fit,gridK)) #On a grid
        standardErrorsGrid_r = robjects.r(''' standardErrorsGrid <-
                                    c(predict.surface.se(fit, extrap=TRUE, nx = '''+str(int(sizePixelMap))+''', ny = '''+str(int(sizePixelMap))+''')) ''')
      else:
        look_r = robjects.r(''' look<-predict.surface(fit, nx = '''+str(int(sizePixelMap))+''', ny = '''+str(int(sizePixelMap))+''') ''')
#standardErrorsGrid <- c(predict.se(fit,gridK)) #On a grid
        standardErrorsGrid_r = robjects.r(''' standardErrorsGrid <-
                                    c(predict.surface.se(fit, nx = '''+str(int(sizePixelMap))+''', ny = '''+str(int(sizePixelMap))+''') ''')
    except: #New version of 'fields' (R package)
      posOutput_surface = 8
      if full:
        look_r = robjects.r(''' look<-predictSurface(fit, extrap=TRUE, nx = '''+str(int(sizePixelMap))+''', ny = '''+str(int(sizePixelMap))+''') ''')
#standardErrorsGrid <- c(predict.se(fit,gridK)) #On a grid
        standardErrorsGrid_r = robjects.r(''' standardErrorsGrid <-
                                    c(predictSurfaceSE(fit, extrap=TRUE, nx = '''+str(int(sizePixelMap))+''', ny = '''+str(int(sizePixelMap))+''')) ''')
      else:
        look_r = robjects.r(''' look<-predictSurface(fit, nx = '''+str(int(sizePixelMap))+''', ny = '''+str(int(sizePixelMap))+''') ''')
#standardErrorsGrid <- c(predict.se(fit,gridK)) #On a grid
        standardErrorsGrid_r = robjects.r(''' standardErrorsGrid <-
                                    c(predictSurfaceSE(fit, nx = '''+str(int(sizePixelMap))+''', ny = '''+str(int(sizePixelMap))+''')) ''')
  gridK_r = robjects.r(''' gridK <- expand.grid(look$x, look$y) ''')
#  linearZ_r = robjects.r(''' linearZ <- expand.grid(look[3]) ''')
  linearZ_r = robjects.r(''' linearZ <- expand.grid(look['''+str(int(posOutput_surface))+''']) ''')
  linearerrZ_r = robjects.r(''' linearerrZ <- expand.grid(standardErrorsGrid['''+str(int(posOutput_surface))+''']) ''')
  tmptab_r = robjects.r(''' tmptab <- cbind(gridK[1],gridK[2],linearZ, linearerrZ) ''')
#
#SAVING in txt
  robjects.r('''filename <- paste(pathOutput,"gridKrig_",label,".txt", sep="")''')
  robjects.r('''write.table(tmptab, filename, sep="\t", col.names = F, row.names = F) ''')
#
  return True


def getAverageDistance(xx, yy, errz=[]): #returns the average distance between the points, weighted by their error
  #
  distances, weights = [], []
  if errz != []:
    for ii in range(len(xx)):
      for jj in range(len(xx)):
        distances.append(numpy.sqrt(((xx[ii]-xx[jj])**2.)+((yy[ii]-yy[jj])**2.)))
        if len(numpy.shape(errz)) == 1: #only one array of errors
          weights.append(1./numpy.sqrt(2.*(errz[ii]**2.)+2.*(errz[jj]**2.)))
        elif len(numpy.shape(errz)) == 2: #only one array of errors
          weights.append(1./numpy.sqrt((errz[0][ii]**2.)+(errz[0][jj]**2.)+(errz[1][ii]**2.)+(errz[1][jj]**2.)))
  #
    return numpy.average(distances, weights=weights)
  else:
    for ii in range(len(xx)):
      for jj in range(len(xx)):
        distances.append(numpy.sqrt(((xx[ii]-xx[jj])**2.)+((yy[ii]-yy[jj])**2.)))
  #
    return numpy.average(distances)

def getMedianDistance(xx, yy): #returns the median distance between the points
  distances = []
  for ii in range(len(xx)):
    for jj in range(len(xx)):
      distances.append(numpy.sqrt(((xx[ii]-xx[jj])**2.)+((yy[ii]-yy[jj])**2.)))
  #
  return numpy.median(distances)


def KrigingMapPython(inputPath, namegal, genTable, label='Z', limits=[-3, +2], visualize=False, sizePixelMap=80):
  #Retrieving galaxy parameters' dictionary
  #Creating the Kriging maps with Python
  #reading input file
  fileKriging = asciidata.open(inputPath+'gridKrig_'+label+'.txt')
  xK, yK, zK, errzK = [], [], [], []
  maxZmap = 0.
  minZmap = 0.
  for jj in range(len(fileKriging[0])):
      xK.append(fileKriging[0][jj])
      yK.append(fileKriging[1][jj])
      if fileKriging[2][jj] != 'NA':
        zK.append(float(fileKriging[2][jj]))
        errzK.append(float(fileKriging[3][jj]))
        if float(fileKriging[2][jj]) > maxZmap: maxZmap = float(fileKriging[2][jj])
        if float(fileKriging[2][jj]) < minZmap: minZmap = float(fileKriging[2][jj])
      else:
        zK.append(nan)
        errzK.append(nan)
  #
  #reshaping
  xK = numpy.array(xK).reshape(sizePixelMap,sizePixelMap)
  yK = numpy.array(yK).reshape(sizePixelMap,sizePixelMap)
  zK = numpy.array(zK).reshape(sizePixelMap,sizePixelMap)
  errzK = numpy.array(errzK).reshape(sizePixelMap,sizePixelMap)

  #
  minZpoints, maxZpoints = numpy.min(genTable[:,2]),  numpy.max(genTable[:,2])
  rangeZmap = [numpy.max([numpy.min([minZpoints, minZmap]), limits[0]]),
               numpy.min([numpy.max([maxZpoints, maxZmap]), limits[1]])]
#
  #
  if savePDF:
    print "Creating Plot"
    if visualize:
      plt.ion()
    else:
      plt.ioff()
    fig = figure(figsize=(6,5))
    clf()
    ax = subplot(111, aspect='equal')
    mapp = ax.pcolor(xK, yK, zK, vmin=rangeZmap[0], vmax = rangeZmap[1])
    ax.set_xlim([numpy.max(xK[0]), numpy.min(xK[0])])
    ax.set_xlabel(r'$\Delta$RA [arcsec]')
    ax.set_ylim([numpy.min(yK), numpy.max(yK)])
    ax.set_ylabel(r'$\Delta$Dec [arcsec]')
#   Isophotes
    from matplotlib.patches import Ellipse
    radiuses = numpy.array([1,3,5,7,9])
    ells = [Ellipse(xy=[0,0], width=(2.*jj*Reff[namegal]/numpy.sqrt(b_a[namegal])),
                height=(2.*jj*Reff[namegal]*numpy.sqrt(b_a[namegal])), angle=90-PA0[namegal],
                edgecolor = 'k', facecolor = 'none', fill = False, linestyle = 'dashed') for jj in radiuses]
    for ee in ells:
      ax.add_artist(ee)
#
#Points
    ax.scatter(numpy.array(genTable[:,0]), numpy.array(genTable[:,1]),
                c=numpy.array(genTable[:,2]),
           vmin=rangeZmap[0], vmax = rangeZmap[1])
  #
    cb = colorbar(mapp)
    if label == 'Z':
      cb.set_label('[Z/H] [dex]')
    elif label == 'CaT':
      cb.set_label(r'CaT index [$\AA$]')
    elif label == 'SN':
      cb.set_label('S/N')
    elif label == 'sigma':
      try:
        cb.set_label(r"$\rm{\sigma}$ [km/s]")
      except:
        cb.set_label("Vel dispersion [km/s]")
    elif label == 'Vel':
      cb.set_label(r"Vel [km/s]")
    ax.set_title(namegal)
  #
    savefig(inputPath+'KrigingMap_python_'+label+'.pdf', bbox_edge = 'tight')
    print "DONE"
  #
  return True


#The code remeasures the Dell with the equation
# angleRot = (numpy.pi/180.)*(90-PA0[ii])
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




def radialProfileLin(namegal, inputFile, label='Z', #binsize=50,  #Bin numerosity
                  binsize=1, #Bin size in arcsec
                  datapoints = []):  #If exist, the radial profiles are limited by the actual datapoints
  #reading input file
  fileKriging = asciidata.open(inputFile)
  xK, yK, zK, errzK = [], [], [], []
  for jj in range(len(fileKriging[0])):
    if fileKriging[2][jj] != 'NA':
      xK.append(fileKriging[0][jj])
      yK.append(fileKriging[1][jj])
      zK.append(float(fileKriging[2][jj]))
      errzK.append(float(fileKriging[3][jj]))
  #
  xK, yK = numpy.array(xK), numpy.array(yK)
  zK, errzK = numpy.array(zK), numpy.array(errzK)
#  xA = -(-xK*numpy.cos((90-PA0[namegal])*numpy.pi/180.) - yK*numpy.sin((90-PA0[namegal])*numpy.pi/180.))
#  yA = -xK*numpy.sin((90-PA0[namegal])*numpy.pi/180.) + yK*numpy.cos((90-PA0[namegal])*numpy.pi/180.)
#  ellDist = numpy.sqrt((b_a[namegal]*(xA**2.))+((yA**2.)/b_a[namegal]))
  ellDist = findDell(xK, yK, PA0[namegal], b_a[namegal])
  ellDist_Sorted = ellDist[permutation_indices(ellDist)]
  zK_Sorted = zK[permutation_indices(ellDist)]
  errzK_Sorted = errzK[permutation_indices(ellDist)]
  #
  # Limit elements within datapoints
  #
  if datapoints != []:
    RA_dp, Dec_dp = numpy.array(datapoints)[:, 0], numpy.array(datapoints)[:, 1]
    ellDist_dp = findDell(RA_dp, Dec_dp, PA0[namegal], b_a[namegal])
    minR, maxR = numpy.min(ellDist_dp), numpy.max(ellDist_dp)
  else:
    minR, maxR = numpy.min(ellDist_Sorted), numpy.max(ellDist_Sorted)
    #
  binR, binZ, bineZ = [], [], []
  for ii in numpy.arange(minR, maxR, binsize):
    tmpR, tmpZ, tmperrZ = [], [], []
    for kk in numpy.arange(len(ellDist_Sorted)):
      if ii <= ellDist_Sorted[kk] < ii+binsize:
        tmpR.append(ellDist_Sorted[kk])
        tmpZ.append(zK_Sorted[kk])
        tmperrZ.append(errzK_Sorted[kk])
    if len(tmpR) > 0:
      binR.append(numpy.average(tmpR))
      binZ.append(numpy.average(tmpZ, weights=1./(numpy.array(tmperrZ)**2.)))
      bineZ.append(numpy.std(tmpZ))
#
  return binR, binZ, bineZ


# In logarithmic bins
def radialProfileLog(namegal, inputFile, label='Z', #binsize=50,  #Bin numerosity
                  binsize=0.01, #Bin size in dex
                  datapoints = []):  #If exist, the radial profiles are limited by the actual datapoints
  #reading input file
  fileKriging = asciidata.open(inputFile)
  xK, yK, zK, errzK = [], [], [], []
  for jj in range(len(fileKriging[0])):
    if fileKriging[2][jj] != 'NA':
      xK.append(fileKriging[0][jj])
      yK.append(fileKriging[1][jj])
      zK.append(float(fileKriging[2][jj]))
      errzK.append(float(fileKriging[3][jj]))
  #
  xK, yK = numpy.array(xK), numpy.array(yK)
  zK, errzK = numpy.array(zK), numpy.array(errzK)
  ellDist = findDell(xK, yK, PA0[namegal], b_a[namegal])
  ellDist_Sorted = ellDist[permutation_indices(ellDist)]
  ellDist_Sorted_log = numpy.log10(ellDist[permutation_indices(ellDist)]/Reff[namegal])
  zK_Sorted = zK[permutation_indices(ellDist)]
  errzK_Sorted = errzK[permutation_indices(ellDist)]
  #
  # Limit elements within datapoints
  #
  if datapoints != []:
    RA_dp, Dec_dp = numpy.array(datapoints)[:, 0], numpy.array(datapoints)[:, 1]
    ellDist_dp_log = numpy.log10(findDell(RA_dp, Dec_dp, PA0[namegal], b_a[namegal])/Reff[namegal])
    minR, maxR = numpy.min(ellDist_dp_log), numpy.max(ellDist_dp_log)
  else:
    minR, maxR = numpy.min(ellDist_Sorted_log), numpy.max(ellDist_Sorted_log)
    #
  binR, binZ, bineZ = [], [], []
  for ii in numpy.arange(minR, maxR, binsize):
    tmpR, tmpZ, tmperrZ = [], [], []
    for kk in numpy.arange(len(ellDist_Sorted_log)):
      if ii <= ellDist_Sorted_log[kk] < ii+binsize:
        tmpR.append(ellDist_Sorted_log[kk])
        tmpZ.append(zK_Sorted[kk])
        tmperrZ.append(errzK_Sorted[kk])
    if len(tmpR) > 0:
      binR.append(numpy.average(tmpR))
      binZ.append(numpy.average(tmpZ, weights=1./(numpy.array(tmperrZ)**2.)))
      bineZ.append(numpy.std(tmpZ))
#
  return binR, binZ, bineZ




def extractStdBin(prof_R, totRealizations, listReal_R, listReal_val):
  prof_errp, prof_errm, prof_median, n_elements = [], [], [], []
  for jj in numpy.arange(len(prof_R)):
    tmpValues = []
    for kk in numpy.arange(totRealizations):
      for ww in numpy.arange(len(listReal_R[kk])):
        #
        if jj == 0: #In case is the first point, all the inner points are included
          if (listReal_R[kk][ww] <= ((prof_R[jj]+prof_R[jj+1])/2.)):
            tmpValues.append(listReal_val[kk][ww])
        #
        elif jj == len(prof_R)-1:#In case is the last point, all the outer points are included
          if (listReal_R[kk][ww] > ((prof_R[jj]+prof_R[jj-1])/2.)):
            tmpValues.append(listReal_val[kk][ww])
        else:
          if (((prof_R[jj]+prof_R[jj-1])/2.) <= listReal_R[kk][ww]
                                      < ((prof_R[jj]+prof_R[jj+1])/2.)):
            tmpValues.append(listReal_val[kk][ww])
    #
    if tmpValues:
      prof_errm.append(numpy.percentile(tmpValues, 16))
      prof_errp.append(numpy.percentile(tmpValues, 84))
      prof_median.append(numpy.percentile(tmpValues, 50))
      n_elements.append(len(tmpValues))
      #
    else:
      prof_errm.append(NaN)
      prof_errp.append(NaN)
      prof_median.append(NaN)
      n_elements.append(0)

  return prof_errp, prof_errm, prof_median, n_elements


####

# For both linear and logarithmic profiles, errors using bootstrapping on the datapoints (if mode = bootstrapping).
# Or, for both the linear and log profiles, errors using MC simulation (if mode = MC)
def MCerrors(linear_prof_R, log_prof_R, totRealizations, namegal, genTable, rangeKriging, label='Z', mode='BS', , sizePixelMap=80.):
  #
  ## create realizations
  list_R_lin, list_Val_lin = [], []
  list_R_log, list_Val_log = [], []
  if not(os.path.exists(pathNick+namegal+'/Kriging/'+mode+'_'+label)):
    os.mkdir(pathNick+namegal+'/Kriging/'+mode+'_'+label)
  #
  for jj in numpy.arange(totRealizations):
    if verbose:
      stdout.write('\r'+mode+'Realization n. %i / %i' % (jj+1, totRealizations))
      stdout.flush()

    if mode == 'BS':
      #
      dummy = bootstrapRealization(genTable, pathNick+namegal+'/Kriging/'+mode+'_'+label, jj)
    elif mode == 'MC':
      dummy = MCextraction(genTable, pathNick+namegal+'/Kriging/'+mode+'_'+label, jj)

    dummy = KrigingR(pathNick+namegal+'/Kriging/'+mode+'_'+label+'/'+mode+str(jj)+'/realization_'+str(jj)+'_Points.txt', visualize=False,
         theta_r = int(rangeKriging), coeff_r = 3, savePdf = False,
         pathOutput = pathNick+namegal+'/Kriging/'+mode+'_'+label+'/'+mode+str(jj)+'/', label=label, sizePixelMap=80.)

  ## extract profiles
  ### LINEAR
    tmpR, tmpVal, tmpeVal = radialProfileLin(namegal, pathNick+namegal+'/Kriging/'+mode+'_'+label+'/'+mode+str(jj)+'/gridKrig_'+label+'.txt',
                                  label=label, datapoints = genTable)
    list_R_lin.append(tmpR)
    list_Val_lin.append(tmpVal)
  #
  ### LOGARITHMIC
    tmpR, tmpVal, tmpeVal = radialProfileLog(namegal, pathNick+namegal+'/Kriging/'+mode+'_'+label+'/'+mode+str(jj)+'/gridKrig_'+label+'.txt',
                                  label=label, datapoints = genTable)
    list_R_log.append(tmpR)
    list_Val_log.append(tmpVal)

  # Extracting standard deviations
  (lin_prof_errp, lin_prof_errm,
        lin_prof_median, n_elements_lin) = extractStdBin(linear_prof_R,
                                                      totRealizations,
                                          list_R_lin, list_Val_lin)
  (log_prof_errp, log_prof_errm,
        log_prof_median, n_elements_log) = extractStdBin(log_prof_R,
                                                      totRealizations,
                                          list_R_log, list_Val_log)

  return (numpy.array(lin_prof_errm), numpy.array(lin_prof_errp),
          numpy.array(lin_prof_median), n_elements_lin,
          array(log_prof_errm), numpy.array(log_prof_errp),
          numpy.array(log_prof_median), n_elements_log)
