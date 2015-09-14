#!/usr/bin/env python

'''
Created on July 14th, 2015

@author: Nicola Pastorello
'''

import glob, os, numpy, time, sys
from sys import stdout
import multiprocessing as mp
from multiprocessing import Process
import platform



######################
# Parallel functions #
######################

def worker(input_queue, output_queue):
  """Start the code processes one after the other"""
  while input_queue.empty() != True:
    listPar = input_queue.get()
    statusReal = job_chunk(listPar)
    if statusReal:
      output_queue.put(True)
    else:
      output_queue.put(False)
  return


def job_chunk(listPar):
    status = runOnGalaxy(listPar)
    return status

def status(proc):
  """Check for processes status"""
  if proc.is_alive==True:
    return 'alive'
  elif proc.is_alive==False:
    return 'dead'
  else:
    return proc.is_alive()


def runOnGalaxy(listPar):
 try:
  import time
  time0 = time.time()
  #Reading old file
  #print ii
  # Create flag file (existing during work in progress)
  open(listPar[1][:-25]+'/inProgress', 'a').close()
  #
  os.system('python '+listPar[1][:59]+'KrigingMapping_v4.py '+listPar[0]+' Z -f')
  #
  #
  print "DONE with "+ii+" in "+str(round((time.time() - time0)/60.,2))+" minutes."
 # os.remove(listPar[1][:-25]+'/inProgress')
  # Create flag file (work done)
 # open(listPar[1][:-25]+'/done', 'a').close()
  #
  return True
#IF ERRORS OCCURED
 except:
  #os.remove(listPar[1][:-25]+'/inProgress')
  #
  return False

########
# MAIN #
########

def mainRun():
  listDirectories = glob.glob('NGC*'); galnames = []
  for ii in listDirectories:
    if not('old' in ii):
      galnames.append(ii)
  dicPathInput = {}
  pathNick = os.getcwd()
  #
  if len(galnames) == 0:
    print "ERROR, NO INPUT DIRECTORIES FOUND"
  else:
    for ii in galnames:
      tmpPath = glob.glob(pathNick+'/'+ii+'/OutputMet_corr_CLEAN.txt')
      dicPathInput[ii] = tmpPath[0]

#Removing galaxies already mapped ('done' flag file in directory)
  for ii in dicPathInput.keys():
    if ((len(glob.glob(ii+'/done')) > 0) or (ii in ['NGC4449', 'NGC5907'])):
      dicPathInput.pop(ii, None)
    if  (len(glob.glob(pathNick+ii+'/inProgress')) > 0):
      os.remove(pathNick+ii+'/inProgress')

  #Let's parallel!
  nproc = mp.cpu_count()
#
  input_queue = mp.Queue()
  output_queue = mp.Queue()
#
  namegals = []
  for ii in dicPathInput.keys():
#  for ii in ['NGC3377']:
    input_queue.put([ii, dicPathInput[ii]])
    namegals.append(ii)
##
  if platform.system() == 'Linux':
    nElements = input_queue.qsize()
  else:
    nElements = len(dicPathInput.keys())
  procs = []    # processes container
  # Start the worker processes
#
  for ii in range(nproc-1):
    print "Process loop ", ii
    procs.append(mp.Process(target=worker, args=(input_queue,output_queue,)))
#
  for ii in procs:
    ii.start()
#
  for ii in procs:
    print "Process ", ii," @ " , ii.pid, " is ", status(ii)
# Wait processes to finish
#    while input_queue.empty() != True:
#      time.sleep(10) # loose 10 seconds
  counter = 0
  while output_queue.qsize() != nElements:
    stdout.write("\rDONE %i/%i %i" % (output_queue.qsize(),
                     nElements, counter))
    stdout.flush()
    counter += 1
    time.sleep(30)
#
  print "FAILED WITH THE GALAXIES:"
  lenOutputQueue = output_queue.qsize()
  if len(namegals) == lenOutputQueue:
    for ii in numpy.arange(lenOutputQueue):
      if not(output_queue.get()):
        print namegals[ii]
  return 'Completed!'



if __name__ == '__main__':
  aa = mainRun()



