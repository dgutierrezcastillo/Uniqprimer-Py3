'''
Created on Jan 1, 2011

@author: John L. Herndon
@contact: herndon@cs.colostate.edu
@organization: Colorado State University
@group: Computer Science Department, Asa Ben-Hur's laboratory 
'''

from . import utils
import tempfile
from . import programs
from . import eprimerparser
from . import primersearchutils
from . import fastaparser





class PrimerManager( object ):
    '''
    A class used to find primers given a set of sequences.
    '''

    def __init__( self, eprimerOptions ):
        self.eprimer = programs.Eprimer( eprimerOptions )
        self.primersearch = programs.PrimerSearch( )
    
                
    def findPrimers( self, sequences, outputFile, primerpairs = 20, returnPrimers = False, chunk_size = 10000 ):
        '''
        A method to find a set of primers based on the given sequences
        '''
        
        utils.logMessage( "PrimerManager::findPrimer(s )", "designing primers for unique sequences" )
        
        # Eliminate all sequences that are less than the minimum product size
        min_prod = 200
        try:
            # Try to get min range from options
            prange = self.eprimer.options.getProductRange()
            if '-' in prange:
                min_prod = int(prange.split('-')[0])
        except:
            pass

        utils.logMessage("PrimerManager::findPrimers()", f"Filtering sequences shorter than {min_prod} bp")
        original_count = len(sequences)
        sequences = [x for x in sequences if len(x) >= min_prod]
        utils.logMessage("PrimerManager::findPrimers()", f"Filtered {original_count - len(sequences)} sequences. {len(sequences)} remain.")

        if not sequences:
            utils.logMessage("PrimerManager::findPrimers()", "No sequences long enough for primer design.")
            if returnPrimers:
                return []
            return

        # Chunk large sequences to avoid eprimer3/primer3 limits
        processed_sequences = []
        for seq in sequences:
            if len(seq) > chunk_size:
                utils.logMessage("PrimerManager::findPrimers()", f"Chunking large sequence of length {len(seq)} with chunk_size {chunk_size}")
                for i in range(0, len(seq), chunk_size):
                    chunk = seq[i : i + chunk_size + 1000] # overlap to avoid missing pairs at boundaries
                    if len(chunk) >= min_prod:
                        processed_sequences.append(chunk)
            else:
                processed_sequences.append(seq)

        utils.logMessage("PrimerManager::findPrimers()", f"Total chunks/sequences to process: {len(processed_sequences)}")

        all_primers = []
        
        # eprimer3 often fails with multi-fasta files or requires special handling.
        # We will process sequences one by one to ensure reliability.
        for i, seq in enumerate(processed_sequences):
            temp_seq_file = utils.getTemporaryDirectory() + f"/seq_{i}.fasta"
            temp_out_file = utils.getTemporaryDirectory() + f"/out_{i}.ep3"
            
            fastaparser.writeFastaFile([seq], temp_seq_file)
            
            try:
                self.eprimer.execute([temp_seq_file, temp_out_file])
                primers = eprimerparser.parsePrimerSequences(temp_out_file)
                all_primers.extend(primers)
                
                # Stop early if we have enough primers to speed things up
                if len(all_primers) >= 100:
                    break
            except Exception as e:
                # utils.logMessage("PrimerManager::findPrimers()", f"Warning: eprimer3 failed for sequence {i}: {e}")
                continue

        utils.logMessage("PrimerManager::findPrimers()", f"Found {len(all_primers)} primers across all sequences.")
        
        if returnPrimers:
            return all_primers
        
    
    def getPrimers( self, sequences, chunk_size=10000 ):
        
        utils.logMessage( "PrimerManager::getCommonPrimers", "finding primers that are common to all include files" )
            
        if len( sequences ) == 0:
            raise utils.NoPrimersExistException( )
        
        referenceEPrimerFile = utils.getTemporaryDirectory( ) + "/referenceprimers.ep3"
        
        #run eprimer to find primers in the reference file
        primers = self.findPrimers( sequences, referenceEPrimerFile, 20, True, chunk_size=chunk_size )
        
        
        if len( primers ) == 0:
             raise utils.NoPrimersExistException( )
        
        return primers
 
    def crossValidatePrimers2( self, primers, includeFile, j ):    
        includeSequences = fastaparser.parseFastaFile( includeFile )
        #write a primer search input file with using the primers argument
        primerInputFileName = utils.getTemporaryDirectory( ) + "/tmpinputprimers2.ps" + str(j)
        primerOutputFileName = utils.getTemporaryDirectory( ) + "/tmpoutputprimers2.ps" + str(j)
        primersearchutils.writePrimerSearchInputFile( primers, primerInputFileName )

        utils.logMessage( "PrimerManager::crossValidatePrimers", "finding primers that are in the supplied include file" )
        #run primer search to identify the primers
        self.primersearch.execute( [ includeFile, primerInputFileName, primerOutputFileName, "0" ] )

        #read the found primers from the file
        commonPrimers = primersearchutils.parsePrimerSearchFile( primerOutputFileName )

        #compose a list of primers that are not found in the exclude file...
        returnValue = [ ]

        for primer in primers:
            if primer.id in commonPrimers:
                returnValue.append( primer )

        utils.logMessage( "PrimerManager::crossValidatePrimers", "{0} unique primers identified out of {1}".format( len( returnValue ), len( primers ) ) )

        if len( returnValue  ) == 0:
            raise utils.NoPrimersExistException( )

        return returnValue

    
    def crossValidatePrimers( self, primers, excludeFile ):             
        
        excludeSequences = fastaparser.parseFastaFile( excludeFile )
        
        #write a primer search input file with using the primers argument
        primerInputFileName = utils.getTemporaryDirectory( ) + "/tmpinputprimers.ps"
        primerOutputFileName = utils.getTemporaryDirectory( ) + "/tmpoutputprimers.ps"
        primersearchutils.writePrimerSearchInputFile( primers, primerInputFileName )
        
        utils.logMessage( "PrimerManager::crossValidatePrimers", "finding primers that are not in the supplied exclude file" )
        #run primer search to identify the primers
        self.primersearch.execute( [ excludeFile, primerInputFileName, primerOutputFileName, "10" ] )
        
        #read the found primers from the file
        commonPrimers = primersearchutils.parsePrimerSearchFile( primerOutputFileName )
        
        #compose a list of primers that are not found in the exclude file...
        returnValue = [ ]
        
        for primer in primers:
            if primer.id not in commonPrimers:
                returnValue.append( primer )
        
        utils.logMessage( "PrimerManager::crossValidatePrimers", "{0} unique primers identified out of {1}".format( len( returnValue ), len( primers ) ) )
        
        if len( returnValue  ) == 0:
            raise utils.NoPrimersExistException( )
        
        return returnValue
    
    
    
