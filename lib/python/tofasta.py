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

import types
import os
import time
import string
import regsub
import regex
import sys

from types import *

# MGI Python libraries

# Sequence Retrieval Tool library

import CGInocontenttype
import EmbossClient
import MouseMineClient
import genomelib
import tofastalib

ec = EmbossClient.EmbossClient(config)

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

    error = ec.error
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
        debug = config.lookup('DEBUG')

        parms = self.get_parms()

        if debug != '0':
            print("Input Parms")
            print(parms)
            print(" ")

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
                
        except (self.error) as message:
            # Give an error screen to the user which passes
            # along the message which was raised with the
            # exception.

            list = [
                '*****',
                'An error occurred while trying to retrieve your ' + \
                'sequence(s) from our EMBOSS repository.',
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
            print(line)

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
    embosssequence = ''
    sequence = ''
    debug = ''
    inputSeqList = []
    seperator = "#SEP#" # 3.4 seperates multilpe entries in one 'seqs' parm

    # clean the input parameters to ensure correct naming of seq parms
    parms = cleanInputParms(parms)

    # check the parameters for error conditions and raise an exception
    # if any are found

    if not ('upfile' in parms or 'seqs' in parms):
        raise ToFASTACGI.error('''Please select at least one sequence.

To make a selection, click the Back button, click the checkbox next to each desired
sequence, and then download in FASTA format.''')

    if 'upfile' in parms and 'seqs' in parms:
        if parms['upfile'] == '' and parms['seqs'] == '':
            raise ToFASTACGI.error('''Input to Sequence Retrieval Tool incorrect.
You have to enter a sequence specification or upload a file.''')
        elif parms['upfile'] != '' and parms['seqs'] != '':
            raise ToFASTACGI.error('''Input to Sequence Retrieval Tool incorrect.
You have entered both a sequence specification and uploaded a file.
Please specify the sequence you wish to retrieve by only one method.''')

    # fill in values from the submitted parameters

    if 'seqs' in parms:
        seqs = parms['seqs']
    if 'returnErrors' in parms:
        returnerrors = string.strip(parms['returnErrors'])
    if 'debug' in parms:
        debug = string.strip(parms['debug'])

    # process seqs to assign values
    if seqs != '':
        if type(seqs) == ListType:

            # test to make sure maximum number of requested sequences not
            # exceeded
            if len(seqs) > string.atoi(config.lookup('MAX_SEQS')):
                raise ToFASTACGI.error('Please contact MGI User Support (mgi-help@informatics.jax.org) to retrieve more than %s sequences.' % config.lookup('MAX_SEQS'))

            # MGI 3.4 release
            # There can now be multiple sequence parameters bound into a 
            # single 'seqs' input field.  This splits them out.
            for inputSeq in seqs:
                for seq in string.split(inputSeq,seperator):
                    inputSeqList.append(seq)

            for seqitem in inputSeqList:
                try:

                # Sept 29, 04
                # Added chromo input parameter
                    [id_db,id,chromo,begin,coorend,strand,flank] = \
                        string.split(seqitem,'!')

                except:
                    raise ToFASTACGI.error('One of the requested sequences %s was not specified properly.  Please resubmit your search.' % seqitem)

                if flank == '':
                    flank = '0'
                if strand == '1':
                    strand = '+'
                elif strand == '0':
                    strand = '-'
                elif strand != '+' and strand != '-' :
                    strand = '+'

                # Sept 29, 04
                # Added to support additional SRT parameter
                if chromo != '' and id_db == 'mousegenome':
                     id = chromo

                # Subtracting 1 from start coordinate.  This is done because
                # of a bug in nigFrag not picking up the first bp
                if begin != '' and begin != '1':
                    origionalBegin = begin
                    newBegin = string.atoi(begin) -1
                    begin = '%s' % newBegin

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

    # process upfile to assign values
    if 'upfile' in parms and upfile == '':
        upfile = parms['upfile']
        outupfile = ''
        upfile_lines = string.split(upfile,'\n')
        # test to make sure maximum number of requested sequences not
        # exceeded
        if len(upfile_lines) > string.atoi(config.lookup('MAX_SEQS')):
            raise ToFASTACGI.error('Please contact MGI User Support (mgi-help@informatics.jax.org) to retreived more than %s sequences.' % config.lookup('MAX_SEQS'))
        for line in upfile_lines:
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
    embossfile,genomeupfile,mousemineFile = separateInputFile(upfile)

    failedgenomemessage = ''
    failedembossmessage = ''

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
                failedseqmsg = "Chr:%s  Begin:%s  End:%s  Flank:%s\n" \
                    "Consider reducing the flanking sequence.\n" \
                    % (chromo,begin,coorend,flank)
                failedgenomemessage = "The Sequence Retrieval Tool failed " + \
                    "to find this genome sequence:\n%s" % failedseqmsg

            # Remove full path from display...  
            if 'GENOMIC_PATH' in config :
                genomesequence = regsub.gsub(config['GENOMIC_PATH'],
                    '',genomesequence) 

    except:
        raise ToFASTACGI.error('Error in retrieving genome build sequences.')

    # Retrieve Non-Genome Build Sequences
    try:
         if embossfile != '':
                         embosssequence, failedembossseqs = ec.getSequences(embossfile)

             # Stamp time
                         profiler.stamp('After call to embosslib.getSequences')

                         if failedembossseqs != '\n':
                                failedembossmessage = "The Sequence Retrieval Tool failed " + \
                                        "to find these sequences:\n%s" % failedembossseqs
                        
    except ec.error as message:
        raise ToFASTACGI.error('Error in retrieving sequences.\n' + message)

    if embosssequence == "":
        sequence = genomesequence
    else:
        sequence = string.rstrip(embosssequence) + "\n" + genomesequence

    # Retrieve strain gene sequences from MouseMine
    mmFasta = ''
    mmFailed = ''
    if mousemineFile:
        ids = []
        for line in mousemineFile.split('\n'):
            fields = line.split('\t')
            if len(fields) > 1:
                ids.append(fields[1])
            
        try:
            mmc = MouseMineClient.MouseMineClient(config['MOUSEMINE_URL'])
            mmFasta, mmFailed = mmc.getFasta(ids)
            if mmFailed:
                mmFailed = mmFailed + '\n'
            
            if sequence:
                sequence = sequence.rstrip() + '\n' + mmFasta
            else:
                sequence = mmFasta

        except:
            mmFailed = 'Retrieval from MouseMine failed (%s): %s' % (ids, str(sys.exc_info()[1]))

    # error reporting
    
    if (failedgenomemessage != '') or (failedembossmessage != '') or (mmFailed != ''):
        print("*****\n" + \
              "An error occurred while trying to retrieve your " + \
              "sequence(s) from our EMBOSS repository.\n-----\n" + \
              failedgenomemessage + failedembossmessage + mmFailed + "*****")

    return sequence,profiler,debug


def separateInputFile(upfile):

    # Purpose: Parse the list of sequences to retrieve and separate
    #       that list in sequences that are from the NCBI genome
    #       build or not, and convert to GCG list files.
    # Returns: Three lists of sequences:
    #    1. sequences to retrieve from NCBI genome build
    #    2. sequences to retrieve from elsewhere (i.e., GCG)
    #    3. sequences to retrieve from MouseMine (strain gene sequences) 
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    gcgupfile = '..\n'
    genomeupfile = ''
    mousemineFile = ''

    for line in string.split(upfile,'\n'): 
        if line != '':
           tokens = string.split(line,'\t')
           if tokens[0] == "mousegenome":
               genomeupfile = genomeupfile + string.joinfields(tokens,'\t') + '\n'
           elif tokens[0] == 'straingene':
               mousemineFile = mousemineFile + string.joinfields(tokens,'\t') + '\n'
           else:
               gcgupfile = gcgupfile + string.joinfields(tokens[:-1],'\t') + '\n'

    # reset 'empty' list files to make testing elsewhere easier
    if gcgupfile == "..\n":
        gcgupfile = ""

    return gcgupfile,genomeupfile,mousemineFile

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

    return id_db

###--------------------------------------------------------------------###
def cleanInputParms(inputParms):
# Purpose: Convert multiple APIs to have same input naming convention
#   Detail page javaScript requires we use more than the 'seqs' parms.
#   seq(n) is now valid, where n is any number (e.g. seq1, seq2, seq3)
#   If a flank(n) parameter is passed, it is appended to coorisponding
#   seq(n)
# Returns: dictionary; like self parms, with seq(n) values now in
#   the seqs parameter 
# Assumes: Nothing
# Effects: Nothing
# Throws:  Nothing

    seqList     = []           # list of strings to become new seqs parm
    upfile      = ''
    flankValues = {}           # temp holding for flank values
    flankValueTemplate = '%s'  # used to cast an int flank value back to 
                               # string, for easy modification
    # regular expressions
    seqReg   = regex.compile('seq[0-9]+')  
    flankReg = regex.compile('flank[0-9]+')

    cgiKeys  = list(inputParms.keys())

    # pull out flanking values from input parms, to be matched to seqN later
    for key in cgiKeys:
        if flankReg.match(key) != -1:
            seqParmNum = key[5:]
            flankValue = flankValueTemplate % (string.atoi(inputParms[key]) * 1000)

            flankValues[seqParmNum] = flankValue

    # gather all values seqN parameter
    for key in cgiKeys:

        # Origional input parameter API spec
        if key == 'seqs':
            if type(inputParms['seqs']) == bytes:
                seqList.append(inputParms['seqs'])
            else:
                for seqsValue in inputParms['seqs']:
                    seqList.append(seqsValue)

        # matched seqN regex
        if seqReg.match(key) != -1:

            # since only string parms can have flanking...
            if type(inputParms[key]) == bytes:

                # Determine if this parameter needs flank appended
                if key[3:] in list(flankValues.keys()):
                    seqList.append(inputParms[key] + flankValues[seqParmNum])
                else:
                    seqList.append(inputParms[key])

        if key == "upfile":
            upfile = inputParms[key]

    newInputParms = {}

    if seqList != []:
        newInputParms['seqs'] = seqList
    if upfile != '':
        newInputParms['upfile'] = upfile

    return newInputParms

