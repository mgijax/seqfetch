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
#    parseParameters(params,profiler)
#    separateInputFile(upfile)
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


from types import *

# MGI Python libraries

import seqlib

# Sequence Retrieval Tool library

import CGInocontenttype
import gcglib
import genomelib
import tofastalib

# Profiler library

import Profiler

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

class ToFASTACGI (CGInocontenttype.CGI):
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

        # Instantiate Profiler class to track time different steps take.
        profiler = Profiler.Profiler()

        # Initialize debug parameter
        debug = '0'

        parms = self.get_parms()

        # Stamp time
        profiler.stamp('After self.get_parms')

        try:
            # The call to parseParameters() may raise the 'error'
            # exception.  We catch it below and display its
            # accompanying message for the user.

            sequence,profiler,debug = parseParameters (parms,profiler)

            # Stamp time
            profiler.stamp('After parseParameters')

            # send the output to the user

            output = [sequence]


        except (self.error, seqlib.error), message:
            # Give an error screen to the user which passes
            # along the message which was raised with the
            # exception.

            list = [
                '*****',
                'An error occurred while trying to retrieve your ' + \
                'sequence(s) from our GCG repository.',
                '-----',
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

        # If in debug mode, then write elapsed time profile.
        if debug == '1':
            profiler.write()

        # Write activity to usage log
        tofastalib.writeToUsageLog(config.lookup('USAGE_LOG'))

        return

###--- Private Functions ---###

def parseParameters (
    parms,        # Dictionary of parameters received from an HTML form,
                  # as returned by CGI.get_parms().
    profiler      # Profiler object to track elapsed time.
    ):
    # Purpose: parse the given set of 'parms' to get and return a tuple
    #    of ten items.  performs error checking to ensure complete and
    #    consistent input.
    # Returns: 1. sequence
    #          2. profilter object
    #          3. debug variable
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
    debug = ''

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
            try:
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

            except:
                raise ToFASTACGI.error, 'Incorrect usage. Be sure to include proper number of fields to identify a sequence and correct number of delimiters.'                

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

    failedgenomemessage = ''
    failedgcgmessage = ''

    # Stamp time
    profiler.stamp('Before call to genomelib.getSequences')

    # Retrieve Genome Build Sequences
    try:
        if genomeupfile != '':
            genomesequence,failedgenomeseqs = \
                genomelib.getSequences (genomeupfile,config)

            # Stamp time
            profiler.stamp('After call to genomelib.getSequences')

            if failedgenomeseqs != '':
                failedgenomemessage = "The Sequence Retrieval Tool failed " + \
                    "to find these genome sequences:\n%s" % failedgenomeseqs

    except:
        raise ToFASTACGI.error, 'Error in retrieving genome build sequences.'


    # Retrieve Non-Genome Build Sequences
    try:
         if gcgupfile != '':
             # Translate the upfile into a properly formatted GCG list file
             gcglistfile = gcglib.generateGCGListFile(gcgupfile)

             # Stamp time
             profiler.stamp('After call to gcglib.generateGCGListFile')

             # Check to see if any sequences have been split by GCG
             checkedlistfile,profiler = \
                 gcglib.checkSeqs(gcglistfile,config,profiler)

             # Stamp time
             profiler.stamp('After call to gcglib.checkSeqs')

             # Run ToFASTA
             gcgsequence,failedgcgseqs,profiler = \
                 gcglib.getSequences (checkedlistfile,gcglistfile,\
                 config,profiler)

             # Stamp time
             profiler.stamp('After call to gcglib.getSequences')

             if failedgcgseqs != '\n':
                failedgcgmessage = "The Sequence Retrieval Tool failed " + \
                    "to find these sequences:\n%s" % failedgcgseqs

    except gcglib.error, message:
        raise ToFASTACGI.error, \
            'Error in retrieving sequences.\n' + message

    if gcgsequence == "":
        sequence = genomesequence
    else:
        sequence = gcgsequence[:-1] + "\n" + genomesequence

    if (failedgenomemessage != '') or (failedgcgmessage != ''):
        print "*****\n" + \
              "An error occurred while trying to retrieve your " + \
              "sequence(s) from our GCG repository.\n-----\n" + \
              failedgenomemessage + failedgcgmessage + "*****"

    return sequence,profiler,debug


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
    elif id_db == "sptrembl":
        id_db = "SWISS-PROT"
    elif id_db == "tigrmgi":
        id_db = "TIGR Mouse Gene Index"
    elif id_db == "tigrrgi":
        id_db = "TIGR Rat Gene Index"
    elif id_db == "tigrhgi":
        id_db = "TIGR Human Gene Index"
    elif id_db == "dotsm":
        id_db = "DoTS Mouse"
    elif id_db == "niamgi":
        id_db = "NIA Mouse Gene Index"
    elif id_db == "mousegenome":
        id_db = "mousegenome"

    return id_db


