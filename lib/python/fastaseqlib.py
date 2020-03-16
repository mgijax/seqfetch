# Name: fastaseqlib.py
# Purpose: provides utility functions for dealing with FASTA formatted
#       sequences
# Assumes:
#       that our PYTHONPATH (sys.path) is set properly for us to find the
#       libraries we import
# Public Functions:
#       appendSequences (firstseq, secondseq, overlap, line_length)
#       reformatSeq (sequence, length, description)
#       getSubSeq (sequence, begin, finish, config)

############
# imports  #
############

# Python libraries

import string
import regsub

####################
# Global Variables #
####################

# exception to be raised in this module
error = "fastaseqlib.py error"

####################
# Public Functions #
####################

def appendSequences (firstseq, secondseq, overlap, line_length):

    # Purpose: Append one FASTA-formatted sequence with another 
    #          allowing for the two input sequences to overlap
    #          by a certain number of base pairs.
    # Returns: A string containing the appended sequence in
    #          FASTA format.
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    # Set overlap length
    overlap_length = string.atoi(overlap)

    #############################################################
    # Step 1:                                                   #
    # Process first sequence:                                   #
    #   1. Separate description line from sequence.             #
    #   2. Clean up sequence.                                   #
    #   3. Remove any portion of sequence that overlaps with    #
    #      second sequence.                                     #
    #############################################################

    # Separate sequence from description line
    first_description = firstseq[0:string.find(firstseq,'\n')+1]
    first_sequenceonly = firstseq[string.find(firstseq,'\n')+1:]

    # Strip out return characters from sequence
    first_sequenceonly = regsub.gsub('\n','',first_sequenceonly)

    # Remove overlapping portion from first_seq 
    first_sequenceonly = first_sequenceonly[:-overlap_length]

    #############################################################
    # Step 2:                                                   #
    # Process second sequence:                                  #
    #   1. Separate description line from sequence.             #
    #   2. Clean up sequence.                                   #
    #############################################################

    # Separate sequence from description line
    second_description = secondseq[0:string.find(secondseq,'\n')+1]
    second_sequenceonly = secondseq[string.find(secondseq,'\n')+1:]

    # Strip out return characters
    second_sequenceonly = regsub.gsub('\n','',second_sequenceonly)

    #############################################################
    # Step 3:                                                   #
    # Append first sequence with second one.                    #
    #############################################################

    # Choose the description line to return depending on whether
    # the first sequence is null or not
    if len(firstseq) == 0:
        description = second_description
    else:
        description = first_description

    # Combine sequence and reformat
    sequence = reformatSeq(first_sequenceonly +\
        second_sequenceonly, line_length, description)

    return sequence

def reformatSeq(sequence, length, description):

    # Purpose: Create a FASTA-formatted sequence given a
    #          sequence, description line, and the number
    #          of sequence letters per line in output.
    # Returns: A string containing the FASTA-format sequence.
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    line_length = string.atoi(length)
    seq_length = len(sequence)

    # Compute number of sequence lines in output
    if line_length > 0:
        num_loops = seq_length/line_length + 1
    else:
        num_loops = 1

    # Construct FASTA-formatted sequence, named new_seq.
    begin = 0
    stop = line_length - 1
    new_seq = description

    for i in range(1, num_loops + 1):
        stop = begin + line_length
        if i == 1:
            new_seq = new_seq + sequence[begin:stop]
        else:
            new_seq = new_seq + "\n" + sequence[begin:stop]
        begin = stop

    return new_seq

def getSubSeq(sequence,begin,finish,config):

    # Purpose: Extract a subsequence from a sequence.
    # Returns: A string containing the subsequence in FASTA format.
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    if (begin != '') and (finish != ''):

        # Convert input coordinates to integers
        coor_begin = string.atoi(begin)
        coor_finish = string.atoi(finish)

        # Separate sequence from description line
        description_line = sequence[0:string.find(sequence,'\n')+1]

        sequence_only = sequence[string.find(sequence,'\n')+1:]

        # Remove all whitespace from sequence
        sequence_only = regsub.gsub('\n','',sequence_only)

        # Get subsequence
        sequence_only = sequence_only[coor_begin - 1:coor_finish]

        # Reformat sequence in FASTA format
        fasta_sequence = reformatSeq(sequence_only, \
            config.lookup('LINE_LENGTH'), description_line)

        sequence = fasta_sequence

    return sequence

