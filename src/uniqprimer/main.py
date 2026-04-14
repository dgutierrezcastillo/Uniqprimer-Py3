#!/usr/bin/python
"""
Uniqprimer - finds primers unique to a genome
"""

import sys
import time
import os
import argparse
from . import utils
from . import includefilemanager
from . import excludefilemanager
from . import primermanager

VERSION = "0.5.4"

class UniqPrimerFinder:
    def __init__(self, include_files, exclude_files, cross_validate, eprimer_options, log_file, output_file, fasta_diff, chunk_size=10000):
        utils.logMessage("UniqPrimerFinder::__init__()", "Initializing UniqPrimerFinder")
        self.include_files = include_files
        self.include_file_manager = includefilemanager.IncludeFileManager()
        
        self.exclude_files = exclude_files
        self.exclude_file_manager = excludefilemanager.ExcludeFileManager()
        
        self.primer_manager = primermanager.PrimerManager(eprimer_options)
        
        self.cross_validate = cross_validate
        self.log_file = log_file
        self.output_file = output_file
        self.fasta_diff = fasta_diff
        self.eprimer_options = eprimer_options
        self.chunk_size = chunk_size

        utils.logMessage("UniqPrimerFinder::__init__()", "Initializing UniqPrimerFinder - complete")
    
    def write_output_file(self, primers, output_file_name, max_results=100):
        """
        Write found primers to output file
        """
        with open(output_file_name, 'w') as f:
            for i, primer in enumerate(primers[:max_results], 1):
                f.write(f"{i}\t{primer.forwardPrimer}\t{primer.reversePrimer}\t{primer.productSize}\n")
        utils.logMessage("UniqPrimerFinder::write_output_file()", f"Output file '{output_file_name}' written.")
            
    def find_primers(self):
        utils.logMessage("UniqPrimerFinder::find_primers()", "Finding primers for include files")
        start_time = time.time()
        
        # Combined exclude files
        utils.printProgressMessage("*** Creating Combined Fasta File for Exclude Files ***")
        for exclude_file in self.exclude_files:
            self.exclude_file_manager.addExcludeFile(exclude_file)
        self.exclude_file_manager.exportSequences()
        
        self.include_file_manager.setExcludeFile(self.exclude_file_manager.getOutputFileName())

        # Determine min_length for sequence extraction based on product size range
        min_prod = 100
        try:
            prange = self.eprimer_options.getProductRange()
            if '-' in prange:
                min_prod = int(prange.split('-')[0])
        except:
            pass

        # Finding unique sequences
        utils.printProgressMessage("*** Finding Sequences Unique to Target Genome ***")
        for include_file in self.include_files:
            self.include_file_manager.processIncludeFile(include_file, min_length=min_prod)
                
        unique_sequences = self.include_file_manager.getUniqueSequences()
        
        utils.printProgressMessage("*** Finding Primers ***")
        # Pass the chunk size to the primer manager
        primers = self.primer_manager.getPrimers(unique_sequences, chunk_size=self.chunk_size)
         
        if self.cross_validate:
            utils.printProgressMessage("*** Cross Validating Primers ***")
            primers = self.primer_manager.crossValidatePrimers(primers, self.exclude_file_manager.getOutputFileName())
            for i, include_file in enumerate(self.include_files, 1):
                primers = self.primer_manager.crossValidatePrimers2(primers, include_file, i)
       		
        utils.logMessage("UniqPrimerFinder::find_primers()", f"found {len(primers)} primers") 
        self.write_output_file(primers, self.output_file)
        
        # Save FASTA for Primer3 if requested
        if self.fasta_diff:
            tmp_fasta = os.path.join(utils.getTemporaryDirectory(), "sequenceForEprimer.fasta")
            if os.path.exists(tmp_fasta):
                import shutil
                shutil.copy2(tmp_fasta, self.fasta_diff)

        end_time = time.time()
        elapsed_min, elapsed_sec = divmod(int(end_time - start_time), 60)
        print(f"*** Time Elapsed: {elapsed_min} minutes, {elapsed_sec} seconds ***")
        print(f"*** Output Written to {self.output_file} ***")

def parse_args():
    parser = argparse.ArgumentParser(description="Uniqprimer - finds primers unique to a genome")
    parser.add_argument("-i", "--include", action="append", required=True, help="Include FASTA file (can be specified multiple times)")
    parser.add_argument("-x", "--exclude", action="append", required=True, help="Exclude FASTA file (can be specified multiple times)")
    parser.add_argument("-o", "--output", default="uPrimer.txt", help="Output file for primers (default: uPrimer.txt)")
    parser.add_argument("-l", "--log", default="log_uniqprimer.txt", help="Log file (default: log_uniqprimer.txt)")
    parser.add_argument("-f", "--fasta", help="Output FASTA of differential sequences")
    parser.add_argument("--productsizerange", default="200-250", help="PCR product size range (default: 200-250)")
    parser.add_argument("--primersize", type=int, default=20, help="Optimal primer size (default: 20)")
    parser.add_argument("--minprimersize", type=int, default=18, help="Minimum primer size (default: 18)")
    parser.add_argument("--maxprimersize", type=int, default=27, help="Maximum primer size (default: 27)")
    
    # New Tm and GC relaxation options
    parser.add_argument("--mintm", type=float, default=57.0, help="Minimum melting temperature (default: 57.0)")
    parser.add_argument("--maxtm", type=float, default=63.0, help="Maximum melting temperature (default: 63.0)")
    parser.add_argument("--opttm", type=float, default=60.0, help="Optimum melting temperature (default: 60.0)")
    parser.add_argument("--mingc", type=float, default=20.0, help="Minimum GC percentage (default: 20.0)")
    parser.add_argument("--maxgc", type=float, default=80.0, help="Maximum GC percentage (default: 80.0)")

    # Customizable chunking
    parser.add_argument("--chunksize", type=int, default=10000, help="Maximum sequence length to send to primer3 at once (default: 10000)")

    parser.add_argument("--crossvalidate", action="store_true", help="Cross validate primers against exclude files")
    parser.add_argument("--keeptempfiles", action="store_true", help="Keep temporary files")
    parser.add_argument("-v", "--version", action="version", version=f"Uniqprimer {VERSION}")

    return parser.parse_args()

def main():
    args = parse_args()
    
    eprimer_options = utils.EPrimerOptions()
    eprimer_options.setProductRange(args.productsizerange)
    eprimer_options.setPrimerSize(args.primersize)
    eprimer_options.setMinPrimerSize(args.minprimersize)
    eprimer_options.setMaxPrimerSize(args.maxprimersize)
    
    # Set Tm and GC options
    eprimer_options.setTm(args.mintm, args.opttm, args.maxtm)
    eprimer_options.setGC(args.mingc, args.maxgc)

    utils.initialize(True, not args.keeptempfiles, args.log)
    
    try:
        utils.logMessage("uniqprimer::main()", f"Include files: {args.include}")
        utils.logMessage("uniqprimer::main()", f"Exclude files: {args.exclude}")
        
        finder = UniqPrimerFinder(
            args.include, args.exclude, args.crossvalidate, 
            eprimer_options, args.log, args.output, args.fasta,
            chunk_size=args.chunksize
        )
        finder.find_primers()
        
    except utils.NoFileFoundException as nfe:
        print(f"Error: File not found: {nfe.filename}")
    except utils.ProgramNotFoundException as pnfe:
        print(f"Error: {pnfe.programName} not found. {pnfe.details}")
    except utils.NoPrimersExistException:
        print("Failure: No unique primers exist for this combination")
        # Ensure output file is created even if empty
        if not os.path.exists(args.output):
            open(args.output, 'a').close()
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
    finally:
        utils.shutdown()
    
    print("*** Finished ***")

if __name__ == '__main__':
    main()
