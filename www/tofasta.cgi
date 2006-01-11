#!./python

# Program: tofasta.cgi
# Purpose: retrieve a set of sequences or subsequences from a particular 
#          set of sequence databases in FASTA format.
# Note: see design document in system_docs directory
# System Requirements Satisfied by This Program:
#   Envvars: none
#    Inputs: Either:
#            1) specification for one sequence
#                1) seqID
#                2) database
#                3) begin coordinate (optional)
#                4) end coordinate (optional)
#            2) specification for a set of sequences
#                1) tab-delimited text file containing these fields
#                   (one line per sequence):
#                    1) database name
#                    2) seqID
#                    3) begin coordinate (optional)
#                    4) end coordinate (optional)
#   Outputs: 1) plain text file of sequences in FASTA format.
#            2) errors specifying when sequences not found in GCG
#               that appear at top of output flanked by '*****'.
#   Exit Codes: none
#   Assumes: 1) Configuration file exists.
#   Other System Requirements:
# Implementation:
#   Modules:
#            tofasta.py - ()
#                Contains routines to retrieve one or more sequences
#                using GCG Wisconsin Package.

############
# imports  #
############

# Python libraries
import sys

if '.' not in sys.path:
    sys.path.insert (0, '.')

# add the MGI standard library directory to the PythonPath so that we can find
# the standard Configuration.py and ignoreDeprecation modules:
MGI_LIBS = '/usr/local/mgi/live/lib/python'
if MGI_LIBS not in sys.path:
	sys.path.insert (0, MGI_LIBS)

# for now, ignore any deprecation errors that could be caused by the
# migration to Python 2.4.2 -- we'll fix them later
import ignoreDeprecation

# MGI Python libraries

#################
# configuration #
#################

import Configuration

# Instantiating a Configuration object will adjust our python path 
# further so that we take our LIBDIRS configuration option into account.
# This helps us find other MGI libraries (including those for ToFASTA).

config = Configuration.get_Configuration ('Configuration', 1)

import tofasta


###########
# main    #
###########

print 'Content-type: text/plain\n'

if __name__ == '__main__':
    mycgi = tofasta.ToFASTACGI()
    mycgi.go()
