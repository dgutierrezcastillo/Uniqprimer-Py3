'''
Created on Jan 1, 2011

@author: John L. Herndon
@contact: herndon@cs.colostate.edu
@organization: Colorado State University
@group: Computer Science Department, Asa Ben-Hur's laboratory 
'''


from . import utils
from . import primersequence

from Bio import SeqIO
from Bio import Seq
from functools import reduce

def parseFastaFileAsPrimerSequence( fileName ):
    
    utils.logMessage("fastaparser::parseFastaFileAsPrimerSequence( )", "parsing fasta file {0}".format( fileName ) )
    returnValue = { }
    
    sequences = SeqIO.parse( open( fileName ), "fasta" )
    
    for sequence in sequences:
        seqdata = primersequence.PrimerSequence( sequence.id, len( sequence ), sequence.seq )
        returnValue[ sequence.id ] = seqdata
    
    utils.logMessage("fastaparser::parseFastaFileAsPrimerSequence( )", "read {0} sequences".format( len( list(returnValue.keys( )) ) ) )
    
    return returnValue
    
def parseFastaFile( fileName ):
    '''
    parse a fasta file and return a list of Bio.Seq
    '''
    utils.logMessage("fastaparser::parseFastaFile( )", "parsing fasta file {0}".format( fileName ) )
    
    sequences =  SeqIO.parse( open( fileName ), "fasta" )
    
    return sequences

def writeFastaFile( sequences, fileName ):
    '''
    write a set of sequences to a fasta file.
    returns the name of the new file
    ''' 
    
    utils.logMessage( "fastaparser::writeFastaFile( )", "Writing {0} sequences to fasta file".format( len( sequences ) ) )
    seqRecords = [ ]
    for i, sequence in enumerate(sequences):
        # Biopython Seq objects or lists of characters need to be converted to strings properly
        if isinstance(sequence, (list, tuple)):
            seqStr = "".join(sequence)
        else:
            seqStr = str(sequence)
            
        seqRecord = SeqIO.SeqRecord( Seq.Seq( seqStr ),  id="seq_{0}".format( i ), description="" )
        seqRecords.append( seqRecord )

    with open( fileName, "w" ) as f:
        SeqIO.write( seqRecords, f, "fasta" )
        
    utils.logMessage( "fastaparser::writeFastaFile( )", "writing fasta file complete" )    
    return fileName
