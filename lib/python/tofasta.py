# Name: tofasta.py
# Purpose: This module represents the real logic behind the tofasta.cgi
#    script.  It is implemented as a module so that we can use a compiled
#    version (tofasta.pyc) rather than needing to compile it each time
#    the CGI script is executed.
# On Import: initializes many global variables
# Assumes: Our PYTHONPATH (sys.path) is set properly so that we can find the
#    Configuration.py module.
# Public Functions:
#    None
# Private Functions:
#    parseParameters(params)
#       separateInputFile(upfile)
# Public Classes:
#    ToFASTACGI
# Sample Usage:
#    mycgi = ToFASTACGI()
#    mycgi.go()

#################
# configuration #
#################

import Configuration
config = Configuration.get_Configuration ('Configuration', 1)

############
# imports  #
############

# Python libraries

import os
import time
import string
import regsub
import regex
import sys

sys.path.insert(0,'/home/blk/public_html/yank/lib/python')
import CGI

from types import *

# MGI Python libraries

import seqlib

# Sequence Retrieval Tool library

import gcglib
import genomelib
import tofastalib

###########################################
# exception values when 'error' is raised #
###########################################

# errors reported to stderr

SEQID_NOT_FOUND = '''An attempt to fetch the sequence with the unique
identifier, %s, from %s has failed.  Please make sure the unique identifier
is correct for the database selected.'''

PARAM_PROBLEM = '''Internal error: This CGI is not passing parameters
correctly to the BLAST program.'''

#############
# classes   #
#############

class ToFASTACGI (CGI.CGI):
    # IS: the ToFASTA CGI script
    # HAS: a set of parameters from an HTML form, including:  id,
    #    id_db.
    # DOES: fetches the sequence from the remote web site and returns
    #       results to the user.

    error = gcglib.error
    genomeerror = genomelib.genomeerror

    def main (self):
        # Purpose: This serves as the (conceptual) main program for
        #    the ToFASTA CGI script.
        # Returns: nothing
        # Assumes: all configuration options are set properly.
        # Effects: may do any or all of the following...  read a
        #    sequence from a sequence database over the web and 
        #    write the result back to the user.
        # Throws: nothing


        parms = self.get_parms()

        try:
            # The call to parseParameters() may raise the 'error'
            # exception.  We catch it below and display its
            # accompanying message for the user.

            sequence = parseParameters (parms)

            # If we got this far, then go ahead and mark up the
            # blast results...

            # send the output to the user

            output = [sequence]


        except (self.error, seqlib.error), message:
            # Give an error screen to the user which passes
            # along the message which was raised with the
            # exception.

            list = [
                '*****',
                'The Sequence Retrieval Tool failed to run.',
                '%s' % message,
                '*****'
                ]

            output = list

            # The call to writeToErrorLog() may raise an IOError
            # which we allow the CGI object's error handler t
            # catch.

            tofastalib.writeToErrorLog (message,config)


        for line in output:
            print line

        # Write activity to usage log
        tofastalib.writeToUsageLog(config.lookup('USAGE_LOG'))

        return

###--- Private Functions ---###

def parseParameters (
    parms        # dictionary of parameters received from an HTML form,
                 # as returned by CGI.get_parms()
    ):
    # Purpose: parse the given set of 'parms' to get and return a tuple
    #    of ten items.  performs error checking to ensure complete and
    #    consistent input.
    # Returns: tuple containing ten items -- (BlastRcd, DatabaseRcd,
    #    sequence, v, b, MaskerRcd, seg, xnu, id, id_db) -- see below
    #    for descriptions of each
    # Assumes: nothing
    # Effects: nothing
    # Throws: 'error' if any problems with the parameters are found
    # Notes:  Parameters expected in 'parms' include...
    #    id - string; seqID or seqID.version to use in retrieving the
    #        actual sequence from the database specified by 'id_db'
    #    id_db - string; identifies the external web-accessible
    #        database from which to retrieve the sequence for 'id'.
    #        The value must be one of those recognized by the
    #        seqlib.getSequence() function.
    #
    
    # set a default for each parameter expected (each described above)

    id = ''
    id_db = ''
    begin = ''
    finish = ''
    strand = '+'
    upfile = ''
    genomesequence = ''
    gcgsequence = ''
    sequence = ''

    # check the parameters for error conditions and raise an exception
    # if any are found

    if not (parms.has_key('upfile') or parms.has_key('seqs')):
        raise ToFASTACGI.error, \
            '''Input to Sequence Retrieval Tool incorrect.
No sequence was specified.'''

    if parms.has_key('upfile') and parms.has_key('seqs'):
        if parms['upfile'] == '' and parms['seqs'] == '':
            raise ToFASTACGI.error, \
                '''Input to Sequence Retrieval Tool incorrect.
You have to enter a sequence specification or upload a file.'''
        elif parms['upfile'] != '' and parms['seqs'] != '':
            raise ToFASTACGI.error, \
                '''Input to Sequence Retrieval Tool incorrect.
You have entered both a sequence specification and uploaded a file.
Please specify the sequence you wish to retrieve by only one method.'''

    # fill in values from the submitted parameters

    if parms.has_key ('seqs'):
        seqs = parms['seqs']
    if parms.has_key ('returnErrors'):
        returnerrors = string.strip(parms['returnErrors'])
    if parms.has_key ('debug'):
	debug = string.strip(parms['debug'])

    # process seqs to assign values
    if seqs != '':
        if type(seqs) == ListType:
            for seqitem in seqs:
                [id_db,id,begin,coorend,strand,flank] = \
                    string.split(seqitem,'!')
                if flank == '':
                    flank = 0
                if strand == '1':
                    strand = '+'
                elif strand == '0':
                    strand = '-'
                else:
                    strand = '+'
                id_db = mapToLogicalDB(id_db)
                if begin != '' and coorend != '':
                    upfile = upfile + "%s\t%s\t%s\t%s\t%s\n" % \
                        (id_db,id,string.atoi(begin)-string.atoi(flank),\
                        string.atoi(coorend)+string.atoi(flank),strand)
                elif begin != '' and coorend == '':
                    upfile = upfile + "%s\t%s\t%s\t\t%s\n" % \
                        (id_db,id,string.atoi(begin)-string.atoi(flank),\
                        strand)
                elif begin == '' and coorend != '':
                    upfile = upfile + "%s\t%s\t\t%s\t%s\n" % \
                        (id_db,id,string.atoi(coorend)+string.atoi(flank),\
                        strand)
                else:
                    upfile = upfile + "%s\t%s\t\t\t%s\n" % (id_db,id,strand)

        else:
            [id_db,id,begin,coorend,strand,flank] = string.split(seqs,'!')
            if flank == '':
                flank = 0
            if strand == '1':
                strand = '+'
            elif strand == '0':
                strand = '-'
            else:
                strand = '+'
            id_db = mapToLogicalDB(id_db)
            if begin != '' and coorend != '':
                upfile = "%s\t%s\t%s\t%s\t%s\n" % \
                    (id_db,id,string.atoi(begin)-string.atoi(flank),\
                    string.atoi(coorend)+string.atoi(flank),strand)
            elif begin != '' and coorend == '':
                upfile = "%s\t%s\t%s\t\t%s\n" % \
                    (id_db,id,string.atoi(begin)-string.atoi(flank),strand)
            elif begin == '' and coorend != '':
                upfile = "%s\t%s\t\t%s\t%s\n" % \
                    (id_db,id,string.atoi(coorend)+string.atoi(flank),strand)
            else:
                upfile = "%s\t%s\t\t\t%s\n" % (id_db,id,strand)

    # process upfile to assign values
    if parms.has_key ('upfile') and upfile == '':
        upfile = parms['upfile']
        outupfile = ''
        for line in string.split(upfile,'\n'):
            if line != '':
                [id_db,id,begin,coorend,strand,flank] = string.split(line,'!')
                id_db = mapToLogicalDB(id_db)
                if begin != '' and coorend != '' and flank != '':
                    outupfile = outupfile + "%s\t%s\t%s\t%s\t%s\n" % \
                        (id_db,id,string.atoi(begin)-string.atoi(flank),\
                        string.atoi(coorend)+string.atoi(flank),strand)
                elif begin != '' and flank != '':
                    outupfile = outupfile + "%s\t%s\t%s\t%s\t%s\n" % \
                        (id_db,id,string.atoi(begin)-string.atoi(flank),\
                        coorend,strand)
                elif coorend != '' and flank != '':
                    outupfile = outupfile + "%s\t%s\t%s\t%s\t%s\n" % \
                        (id_db,id,begin,\
                        string.atoi(coorend)+string.atoi(flank),strand)
                else:
                    outupfile = outupfile + "%s\t%s\t%s\t%s\t%s\n" % \
                        (id_db,id,begin,coorend,strand)
        upfile = outupfile

    # Separate upfile into lists of non-genome-build sequences
    # and genome-build sequences
    gcgupfile,genomeupfile = separateInputFile(upfile)

    failedmessage = ''

    # Retrieve Genome Build Sequences
    try:
        if genomeupfile != '':
            genomesequence,failedgenomeseqs = \
                genomelib.getSequences (genomeupfile,config)

            if failedgenomeseqs != '':
                failedmessage = "Sequence Retrieval Tool Failed " + \
                    "to find these genome sequences:\n%s" % failedgenomeseqs

    except:
        raise ToFASTACGI.error, 'Error in retrieving genome build sequences.'


    # Retrieve Non-Genome Build Sequences
    try:
         if gcgupfile != '':
             # Translate the upfile into a properly formatted GCG list file
             gcglistfile = gcglib.generateGCGListFile(gcgupfile)

             # Check to see if any sequences have been split by GCG
             checkedlistfile = gcglib.checkSeqs(gcglistfile,config)

             # Run ToFASTA
             gcgsequence,failedgcgseqs = \
                 gcglib.getSequences (checkedlistfile,gcglistfile,config)

             if failedgcgseqs != '\n':
                failedmessage = failedmessage + \
                    "Sequence Retrieval Tool Failed to find these " +\
                    "non-genome build sequences:\n%s" % failedgcgseqs

    except gcglib.error, message:
        raise ToFASTACGI.error, \
            'Error in retrieving non-genome build sequences.\n' + message

    sequence = gcgsequence[:-1] + "\n" + genomesequence

    if failedmessage != '':
        print "*****\n" + failedmessage + "*****\n"

    return sequence


def separateInputFile(upfile):

    # Purpose: Parse the list of sequences to retrieve and separate
    #       that list in sequences that are from the NCBI genome
    #       build or not, and convert to GCG list files.
    # Returns: Two lists of sequences:
    #       1. sequences to retrieve from NCBI genome build
    #          2. sequences to retrieve from elsewhere (i.e., GCG)
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    gcgupfile = '..\n'
    genomeupfile = ''

    for line in string.split(upfile,'\n'): 
        if line != '':
           tokens = string.split(line,'\t')
           if tokens[0] == "mousegenome":
               genomeupfile = genomeupfile + \
                   string.joinfields(tokens,'\t') + '\n'
           else:
               gcgupfile = gcgupfile + \
                   string.joinfields(tokens[:-1],'\t') + '\n'

    # reset 'empty' list files to make testing elsewhere easier
    if gcgupfile == "..\n":
        gcgupfile = ""

    return gcgupfile,genomeupfile

def mapToLogicalDB(id_db):

    if id_db == "genbank":
        id_db = "Sequence DB"
    elif id_db == "refseq":
        id_db = "RefSeq"
    elif id_db == "swissprot":
        id_db = "SWISS-PROT"
    elif id_db == "trembl":
        id_db = "TrEMBL"
    elif id_db == "tigrmgi":
        id_db = "TIGR Mouse Gene Index"
    elif id_db == "tigrrgi":
        id_db = "TIGR Rat Gene Index"
    elif id_db == "tigrhgi":
        id_db = "TIGR Human Gene Index"
    elif id_db == "dots":
        id_db = "DoTS"
    elif id_db == "mousegenome":
        id_db = "mousegenome"

    return id_db




