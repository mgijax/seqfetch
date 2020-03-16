# Name: tofastalib.py
# Purpose: provides utility functions for the ToFASTA program
# Assumes:
#       that our PYTHONPATH (sys.path) is set properly for us to find the
#       libraries we import
# Public Functions:
#       writeToErrorLog(s,config)
#       writeToUsageLog (
#        program,        # string; BLAST program to run
#        database,       # string; sequence database to search
#        masker,         # string; type of masking to do
#        length,         # integer; length of the sequence
#        id,             # string; the sequence ID
#        id_db,          # string; ID of the sequence database type
#        num_results,    # integer; number of results found
#        num_sequences,  # integer; number of sequences in the database
#        db_size,        # integer; number of letters in the sequence database
#        elapsed         # float; elapsed time to run BLAST
#        ):

############
# imports  #
############

# Python libraries

import os
import time

####################
# Global Variables #
####################

# exception to be raised in this module
error = "tofastalib.py error"

####################
# Public Functions #
####################

def writeToErrorLog (
        s,               # string; line to write to the error log
        config
        ):
        # Purpose: write 's' to MouseBLAST's error log, along with the IP
        #       address of the client and the current time
        # Returns: nothing
        # Assumes: the error log is writeable
        # Effects: appends to the file at 'ERRORLOG'
        # Throws: IOError if we cannot write to the file

        ip = 'unknown'
        if 'REMOTE_ADDR' in os.environ:
                ip = os.environ['REMOTE_ADDR']
        timestamp = time.asctime (time.localtime(time.time()))

        entry = '\n%s\t%s\t%s' % (timestamp, ip, s)

        fd = open (config.lookup('ERROR_LOG'), 'a')
        fd.write (entry)
        fd.close ()
        return

def writeToUsageLog (logfilename):
#        program,        # string; BLAST program to run
#        database,       # string; sequence database to search
#        masker,         # string; type of masking to do
#        length,         # integer; length of the sequence
#        id,             # string; the sequence ID
#        id_db,          # string; ID of the sequence database type
#        num_results,    # integer; number of results found
#        num_sequences,  # integer; number of sequences in the database
#        db_size,        # integer; number of letters in the sequence database
#        elapsed         # float; elapsed time to run BLAST
#        ):
        # Purpose: write a line to MouseBLAST's usage log, containing the
        #       parameters as well as the client's IP address and the current
        #       time
        # Returns: nothing
        # Assumes: the file at 'USAGELOG' is writeable
        # Effects: appends a line to MouseBLAST's usage log
        # Throws: IOError if we cannot write to the file

        ip = 'unknown'
        if 'REMOTE_ADDR' in os.environ:
                ip = os.environ['REMOTE_ADDR']

        timestamp = time.asctime (time.localtime(time.time()))

        entry = '\n%s\t%s' % \
                (timestamp, ip)

        fd = open (logfilename, 'a')
        fd.write (entry)
        fd.close ()
        return

def writeToErrorLogDebug(s,     config):
        if config.lookup('DEBUG') == '1':
                writeToErrorLog ('DEBUG: ' + s, config)
