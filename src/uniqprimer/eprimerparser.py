'''
Created on Jan 1, 2011

@author: John L. Herndon
@contact: herndon@cs.colostate.edu
@organization: Colorado State University
@group: Computer Science Department, Asa Ben-Hur's laboratory 
'''

import os
from . import utils
import re

def parsePrimerSequences( p3File ):
    '''
    parse a primer3 boulder output file for all primers
    '''
    
    utils.logMessage( "eprimerparser::parsePrimerSequences( )", "parsing for primer sequences" )
    if os.path.exists( p3File ) == False:
        utils.logMessage( "eprimerparser::parsePrimerSequences( )", "ERROR - primer3 output file was not found" )
        raise utils.NoFileFoundException( p3File )
    
    results = {}
    with open(p3File, 'r') as f:
        for line in f:
            if '=' in line:
                key, val = line.strip().split('=', 1)
                results[key] = val
    
    num_returned = int(results.get('PRIMER_PAIR_NUM_RETURNED', 0))
    primers = []
    
    for i in range(num_returned):
        primer_set = utils.PrimerSet(str(i))
        
        # Extract data from primer3 output keys
        fwd_seq = results.get(f'PRIMER_LEFT_{i}_SEQUENCE')
        rev_seq = results.get(f'PRIMER_RIGHT_{i}_SEQUENCE')
        prod_size = results.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE')
        fwd_tm = results.get(f'PRIMER_LEFT_{i}_TM')
        rev_tm = results.get(f'PRIMER_RIGHT_{i}_TM')
        
        if fwd_seq and rev_seq:
            primer_set.setForwardPrimerData(fwd_seq, fwd_tm)
            primer_set.setReversePrimerData(rev_seq, rev_tm)
            if prod_size:
                primer_set.setProductSize(int(prod_size))
            primers.append(primer_set)
    
    utils.logMessage( "eprimerparser::parsePrimerSequences( )", "finished parsing. found {0} primers".format( len( primers ) ) )
    return primers        
