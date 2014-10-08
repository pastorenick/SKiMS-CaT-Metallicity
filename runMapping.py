#!/usr/bin/env python

'''
Created on April 11th, 2014

@author: Nicola Pastorello
'''
# Without sigma
#from KrigingMapping_parallel import *
# With sigma
from KrigingMapping_parallel_v2 import *
 
__builtins__.thetaFromDic = True #Instead of measuring the average distance between the points, 
                    # it takes the theta from the dictionary
__builtins__.pathNick = '/nfs/cluster/gals/npastore/exLustre/REDUCTION/SKiMS_CaTmet/'
__builtins__.savePDF = False	#In queue it doesn't work
dummy = runMapping()

