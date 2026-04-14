'''
Created on Jan 1, 2011

@author: John L. Herndon
@contact: herndon@cs.colostate.edu
@organization: Colorado State University
@group: Computer Science Department, Asa Ben-Hur's laboratory 
'''

from . import utils
import os
import subprocess


class ProgramBase:
    def __init__(self):
        self.programName = None
    
    def getProcessArgs(self, args):
        raise NotImplementedError("Subclasses must implement getProcessArgs")
        
    def execute(self, args, is_async=False):
        """
        Run an external program
        """
        utils.logMessage("ProgramBase::execute()", f"Running the {self.programName} program.")
        args, outputFile = self.getProcessArgs(args)
        utils.printProgressMessage(f"*** Running {self.programName} ***")
        utils.logList("ProgramBase::execute()", args)
        
        try:
            if is_async:
                subprocess.Popen(args)
                return outputFile
            
            # For primer3_core, we need to handle stdin
            if self.programName == "primer3_core":
                # The first arg in inputArgs is actually the input string for primer3
                input_str = args[0]
                result = subprocess.run([self.primer3core], input=input_str, capture_output=True, text=True, check=True)
                with open(outputFile, 'w') as f:
                    f.write(result.stdout)
            else:
                result = subprocess.run(args, capture_output=True, text=True, check=True)
            
            utils.logMessage("ProgramBase::execute()", f"{self.programName} output captured to {outputFile}")
            if result.stderr:
                utils.logMessage("ProgramBase::execute()", f"{self.programName} errors: {result.stderr}")
            
            utils.printProgressMessage(f"*** Running {self.programName} Complete ***")
            return outputFile
            
        except subprocess.CalledProcessError as e:
            utils.logMessage("ProgramBase::execute()", f"Error running {self.programName}: {e.stderr}")
            raise Exception(f"{self.programName} failed with return code {e.returncode}: {e.stderr}")
        except Exception as e:
            utils.logMessage("ProgramBase::execute()", f"Unexpected error running {self.programName}: {e}")
            raise

class Nucmer(ProgramBase):
    def __init__(self):
        super().__init__()
        nucmerPath = utils.search_file('nucmer')
        if nucmerPath is None:
            raise utils.ProgramNotFoundException('nucmer', "Please ensure that the MUMmer package is installed.")
        
        self.nucmer = nucmerPath
        self.programName = "nucmer"
        
    def getProcessArgs(self, inputArgs):
        identifier = "nucmer_alignments"
        args = [self.nucmer, '-p', identifier, '--coords', '--minmatch', '300', '--maxgap', '1']
        args.extend(inputArgs)
        outputFile = f"{identifier}.coords"
        return args, outputFile

class Eprimer(ProgramBase):
    """
    Refactored to call primer3_core directly instead of eprimer3 wrapper for better control.
    """
    def __init__(self, eprimerOptions):
        super().__init__()
        self.programName = "primer3_core"
        self.options = eprimerOptions
        
        primer3corePath = utils.search_file("primer3_core")
        if primer3corePath is None:
            raise utils.ProgramNotFoundException("primer3_core", "Please ensure that the primer3 package is installed.")
        self.primer3core = primer3corePath
        
    def getProcessArgs(self, inputArgs):
        # inputArgs[0] = fasta file containing one sequence
        # inputArgs[1] = output file name
        
        fastaFile = inputArgs[0]
        outputFile = inputArgs[1]
        
        # Read sequence from FASTA
        sequence = ""
        with open(fastaFile, 'r') as f:
            for line in f:
                if not line.startswith(">"):
                    sequence += line.strip()
        
        # Build primer3 boulder format input
        p3_input = [
            f"SEQUENCE_ID=seq",
            f"SEQUENCE_TEMPLATE={sequence}",
            "PRIMER_TASK=generic",
            "PRIMER_PICK_LEFT_PRIMER=1",
            "PRIMER_PICK_INTERNAL_OLIGO=0",
            "PRIMER_PICK_RIGHT_PRIMER=1",
            f"PRIMER_OPT_SIZE={self.options.primerSize}",
            f"PRIMER_MIN_SIZE={self.options.minPrimerSize}",
            f"PRIMER_MAX_SIZE={self.options.maxPrimerSize}",
            f"PRIMER_PRODUCT_SIZE_RANGE={self.options.productRange}",
            f"PRIMER_MIN_TM={self.options.minTm}",
            f"PRIMER_OPT_TM={self.options.optTm}",
            f"PRIMER_MAX_TM={self.options.maxTm}",
            f"PRIMER_MIN_GC={self.options.minGC}",
            f"PRIMER_MAX_GC={self.options.maxGC}",
            "PRIMER_NUM_RETURN=5",
            "="
        ]
        
        return ["\n".join(p3_input)], outputFile
    
class PrimerSearch(ProgramBase):
    def __init__(self):
        super().__init__()
        self.programName = "PrimerSearch"
        primerSearchPath = utils.search_file("primersearch")
        if primerSearchPath is None:
            raise utils.ProgramNotFoundException("primersearch", "Please ensure that the EMBOSS package is installed.")
        self.primerSearch = primerSearchPath
        
    def getProcessArgs(self, inputArgs):
        args = [
            self.primerSearch,
            '-seqall', inputArgs[0],
            '-infile', inputArgs[1],
            '-mismatchpercent', inputArgs[3],
            '-outfile', inputArgs[2]
        ]
        return args, inputArgs[2]
