'''
Created on Jan 1, 2011

@author: John L. Herndon
@contact: herndon@cs.colostate.edu
@organization: Colorado State University
@group: Computer Science Department, Asa Ben-Hur's laboratory 
'''

import time
import os.path
from os import pathsep
import tempfile
import shutil
from functools import reduce

def getTimeStamp():
    return time.strftime('%d%m%Y-%H%M%S')


class Match:
    '''
    Record where two genomes line up. Stores only alignments for one part of the genome.
    '''
    def __init__(self, start, end, seqID):
        self.seqID = seqID
        self.start = start
        self.end = end
    
    def __repr__(self):
        return f"Start: {self.start}, End: {self.end}, SeqID: {self.seqID}"


class PrimerSet:
    def __init__(self, id):
        self.id = id
        self.productSize = 0
        self.forwardPrimer = ""
        self.forwardMeltTemp = ""
        self.reversePrimer = ""
        self.reverseMeltTemp = ""
        
    def setProductSize(self, productSize):
        self.productSize = productSize
        
    def setForwardPrimerData(self, sequence, temp):
        self.forwardPrimer = sequence
        self.forwardMeltTemp = temp
        
    def setReversePrimerData(self, sequence, temp):
        self.reversePrimer = sequence
        self.reverseMeltTemp = temp
    

# Search function from http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
def search_file(filename):
    """ Find file """
    search_path = os.getenv('PATH')
    logMessage("utils::search_file", f"Path: {search_path}")
    file_found = False
    paths = search_path.split(pathsep)
    for path in paths:
        if os.path.exists(os.path.join(path, filename)):
            file_found = True
            break
    if file_found:
        return os.path.abspath(os.path.join(path, filename))
    else:
        return None

tempDir = ""
removeTemp = True
verbose = False
logFile = None

def initialize(isVerbose, cleanup, lf):
    global removeTemp
    global tempDir
    global verbose
    global logFile

    logFile = lf
    verbose = isVerbose
    tempDir = tempfile.mkdtemp()
    initializeLogging()
    removeTemp = cleanup
    logMessage("utils::Initialize()", f"Initialization complete. Temporary directory: {tempDir}")
    
def printProgressMessage(message):
    global verbose
    if verbose:
        print(message)

def getTemporaryDirectory():
    global tempDir
    return tempDir

def initializeLogging():
    global logFile
    if logFile:
        logFile = open(logFile, 'w')
    
def shutdown():
    global removeTemp
    global tempDir
    shutdownLogging()
    if removeTemp and tempDir and os.path.exists(tempDir):
        print("*** Removing temporary directory ***")
        shutil.rmtree(tempDir)
    
def shutdownLogging():
    global logFile
    if logFile is not None:
        logFile.close()
        logFile = None

def logList(method, lst):
    message = reduce(lambda x, y: f"{x} {y}", lst)
    logMessage(method, message)
    
def logMessage(method, message):
    global logFile
    if logFile is None:
        return
    log = f"{method} - {message}"
    logFile.write(log + "\n")
    logFile.flush()

class EPrimerOptions:
    def __init__(self):
        self.minPrimerSize = 18
        self.maxPrimerSize = 27
        self.primerSize = 20
        self.productRange = "200-250"
        self.minTm = 57.0
        self.maxTm = 63.0
        self.optTm = 60.0
        self.minGC = 20.0
        self.maxGC = 80.0
    
    def setPrimerSize(self, size):
        size = int(size)
        if size > 35:
            size = 35
        
        self.primerSize = size
        if self.primerSize < self.minPrimerSize:
            self.maxPrimerSize = self.primerSize 
        elif self.primerSize > self.maxPrimerSize:
            self.maxPrimerSize = self.primerSize
    
    def getPrimerSize(self):
        return self.primerSize
    
    def setMinPrimerSize(self, minSize):
        self.minPrimerSize = minSize
    
    def getMinPrimerSize(self):
        return self.minPrimerSize

    def setMaxPrimerSize(self, size):
        self.maxPrimerSize = size
        
    def getMaxPrimerSize(self):
        return self.maxPrimerSize
    
    def setProductRange(self, range):
        self.productRange = range
        
    def getProductRange(self):
        return self.productRange

    def setTm(self, min_tm, opt_tm, max_tm):
        self.minTm = float(min_tm)
        self.optTm = float(opt_tm)
        self.maxTm = float(max_tm)

    def setGC(self, min_gc, max_gc):
        self.minGC = float(min_gc)
        self.maxGC = float(max_gc)

class NoPrimersExistException(Exception):
    def __init__(self):
        super().__init__()

class ProgramNotFoundException(Exception):
    def __init__(self, programName, details):
        super().__init__()
        self.programName = programName
        self.details = details
        
class NoFileFoundException(Exception):
    def __init__(self, filename):
        super().__init__()
        self.filename = filename
        
class ModuleNotInitializedException(Exception):
    def __init__(self, moduleName, reason):
        super().__init__()
        self.moduleName = moduleName
        self.reason = reason
