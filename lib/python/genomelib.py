# Name: genomelib.py
# Purpose: provides utility functions for dealing with NCBI genome builds.
# Assumes:
#       that our PYTHONPATH (sys.path) is set properly for us to find the
#       libraries we import
# Public Functions:
#       getSequences (listfile, config)
# Private Functions:
#       runNames (listfile, filename, gcgserver, login)
#       runToFASTA2 (cbrout, filename, gcgserver, login)

############
# imports  #
############

# Python libraries

import tempfile
import string
import os

# MGI Python libraries

import runCommand

####################
# Global Variables #
####################

# exception to be raised in this module
genomeerror = "genomelib.py error"

####################
# Public Functions #
####################

def getSequences (
        upfile,
        config                  # configuration object
        ):

    # Purpose: Retrieves a subsequence from the mouse genome build.
    #          The chromosome, begin and end coordinates, and the
    #          strand are required.
    # Returns: The genome subsequence in FASTA format with a description
    #          line containing information about the chromosome and
    #          subsequence coordinates.
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    try:

        sequence = ''
        sequences = ''
        seqs = []
        failedseqs = ''

        # Process list file containing specifications for subsequences
        # to retrieve.

        for line in string.split(upfile,'\n'):
            # chromsome number
            chromosome = ''

            # subsequence begin coordinate	
            begin = ''

            # subsequence end coordinate
            ends = ''

            # subsequence strand ('+' or '-')
            strand = ''

            if line != "":        
                tokens = string.split(line,'\t')
                if len(tokens) == 5:
                    chromosome = tokens[1]
                    if tokens[2]:
                        if tokens[2] != '':
                            begin = string.atoi(tokens[2])
                    if tokens[3]:
                        if tokens[3] != '':
                            ends = string.atoi(tokens[3])
                    if tokens[4]:
                        if tokens[4] != '':
                            strand = tokens[4]

                seqs.append([chromosome,begin,ends,strand])

        # Retrieve subsequences

        for seq in seqs:

            chromosome = ''
            begin = ''
            stop = ''
            strand = '+'

            if seq[0]:
                chromosome = seq[0]
            if seq[1]:
                begin = seq[1]
            if seq[2]:
                stop = seq[2]
            if seq[3]:
                strand = seq[3]

            outputfile = tempfile.mktemp()

            cmd = config.lookup('NIBFRAGDIR') + "/nibFrag " +\
                config.lookup('NIBDBS') +\
                "/%s.fa.nib %s %s %s %s" %\
                (chromosome,begin,stop,strand,outputfile)

            stdout,stderr,exitcode = runCommand.runCommand(cmd, {})

            if exitcode != 0:
                failedseqs = failedseqs + \
                    "mousegenome!%s!%s!%s!%s!\n" % \
                    (chromosome,begin,stop,strand)
            else:
                outfd = open(outputfile,'r')
                outfilelines = outfd.readlines()
                outfd.close()

                # remove file
                os.unlink(outputfile)

                sequences = sequences + string.joinfields(outfilelines,'')

    except genomeerror, message:
        raise genomeerror, '' + message

    return sequences, failedseqs
