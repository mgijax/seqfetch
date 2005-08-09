# Name: gcglib.py
# Purpose: provides utility functions for dealing with GCG programs that
#       are run on a remote server.
# Assumes:
#	that our PYTHONPATH (sys.path) is set properly for us to find the
#	libraries we import
# Public Functions:
#	generateGCGListFile (tablist)
#	checkSeqs (input, config, profiler)
#	getSequences (listfile, firstlistfile, config, profiler)
# Private Functions:
#	runNames (listfile, filename, gcgserver, login)
#	runToFASTA (cbrout, filename, gcgserver, login)
#	getGCGOutput (filename, gcgserver, login, password)
#	copyListFile (filename, gcgserver, login, password)
#	seqReader (loopseq, config)
#	mergeListFiles (nameslistfile, listfile)
#	writeFile (content, config)
#       processSequences (sequence, seqids, splitsubseqs, config):

############
# imports  #
############

# Python libraries

import tempfile
import string
import regex
import regsub
import ftplib
import os

# MGI Python libraries

import runCommand
import seqlib
import fastaseqlib
import tofastalib


####################
# Global Variables #
####################

# exception to be raised in this module
error = "gcglib.py error"

# String that is used by the private printString function to save
# the contents of a file on the remote GCG server using ftplib.
# The variable is reset before each file is retrieved.
contents = ''


####################
# Public Functions #
####################

def getSequences (
    filename,        # string; the filename of GCG list file on GCG server
    firstlistfile,   # string: contents of original list file not created
                     #    by GCG Names command.
    config,          # configuration object
    profiler         # Profiler object to track elapsed time.
    ):

    # Purpose: Retrieve a set of sequences in FASTA-format specified 
    #          in a GCG-format list file.  The list file specifies 
    #          the sequence database, identifier and optionally the
    #          coordinates of a subsequence.
    # Returns: 1. A string containing FASTA-format sequences.
    #          2. List of sequences that were not retrieved.
    #          3. Profiler object for tracking elapsed time.
    # Assumes: 1. Only one copy of a sequence exists in any given
    #             sequence database except for GenBank.  This implies
    #             that all sequence databases in GCG are maintained by
    #             a non-incremental update process (i.e. delete
    #             and reload) except for GenBank.
    #          2. At most two copies of a GenBank sequence exists
    #             in GCG.  One copy from the release, and one from
    #             the update process (GCG logical 'gb_new').
    # Effects: nothing
    # Throws:  nothing

    try:

        #############################################################
        # Step 1:                                                   #
        # Merge list file from GCG Names with original list file to #
        # generate a new list file that:                            #
        #   1. specifies the newest copy of a GenBank sequence if   #
        #      there are more than one copy in the GenBank GCG      #
        #      database (even for split sequences).                 #
        #   2. has subsequence coordinates if sequences not split   #
        #      since these are stripped away while running GCG      #
        #      Names.                                               #
        #   3. lists all split sequences (without subsequence       #
        #      coordinates).                                        #
        #############################################################

        #############################################################
        # Step 1.1. Get contents of output file from GCG Names on   #
        #           remote GCG server.                              #
        #############################################################

        # Stamp time
        profiler.stamp('Before call to getGCGOutput to get output' +\
            ' from GCG Names')

        listfile = getGCGOutput(filename,config.lookup('GCG_SERVER'),
            config.lookup('GCG_ACCOUNT'),config.lookup('GCG_PASS'))

        # Stamp time
        profiler.stamp('After call to getGCGOutput to get output' +\
            ' from GCG Names')

        #############################################################
        # Step 1.2. Merge original GCG list file with GCG list file #
        #           from GCG names as described above.              #
        #############################################################

        mergedlistfile, seqids, splitsubseqs, failedseqs = \
            mergeListFiles(listfile, firstlistfile)

        # Stamp time
        profiler.stamp('After call to mergeListFiles')

        if mergedlistfile == "..\n":
            raise error, 'All sequences missing'

        #############################################################
        # Step 1.3. Write file containing merged list file and copy #
        #           to GCG server.                                  #
        #############################################################

        # Create a list file containing the contents of the 'listfile'
        # variable
        mergedfilename = writeFile(mergedlistfile, config)

        # Copy list file to GCG server
        copyListFile (mergedfilename,config.lookup('GCG_SERVER'),\
            config.lookup('GCG_ACCOUNT'),config.lookup('GCG_PASS'))

        # Stamp time
        profiler.stamp('After call copyListFile to put list file' +\
            ' on GCG server')

        # Remove list file from local server
        os.unlink(mergedfilename)

        #############################################################
        # Step 2:                                                   #
        # Get sequences in merged list file using GCG ToFASTA.      #
        #############################################################

        # Create temporary filename for output of ToFASTA
        tempfile.tempdir = config.lookup('GCG_TEMP_DIR')
        cbrout = tempfile.mktemp()

        # Stamp time
        profiler.stamp('Before call to runToFASTA')

        # Actually run ToFASTA
        runToFASTA (cbrout,mergedfilename,config.lookup('GCG_SERVER'),\
            config.lookup('GCG_ACCOUNT'))

        # Stamp time
        profiler.stamp('After call to runToFASTA')

        # Get sequences from ToFASTA
        sequences = getGCGOutput (cbrout, \
            config.lookup('GCG_SERVER'),config.lookup('GCG_ACCOUNT'),\
            config.lookup('GCG_PASS'))

        # Stamp time
        profiler.stamp('After call to getGCGOutput')


        #############################################################
        # Step 3:                                                   #
        # Process sequences from GCG ToFASTA to:                    #
        #     1. return one sequence per seqID and check for cases  #
        #        where sequence not returned.                       #
        #     2. join together any split sequences.                 #
        #     3. get subsequence of any split sequences if          #
        #        specified in original list file.                   #
        #############################################################

        sequences = processSequences (sequences, seqids, splitsubseqs,
            config)

        profiler.stamp('After call to processSequences')

    except error, message:
        raise error, '' + message

    return sequences, failedseqs, profiler

def checkSeqs(input,config,profiler):

    # Purpose: Runs the program GCG Names to determine whether any
    #          input sequences have been split up by GCG.  Sequences
    #          longer than 350kb are split up into 110kb fragments
    #          that overlap by 10kb.
    # Returns: 1. The filename of the output from GCG Names on the 
    #             remote GCG server.
    #          2. Profiler object for tracking elapsed time.
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    #############################################################
    # Step 1:                                                   #
    # Construct input list file for GCG Names.                  #
    #############################################################

    new_listfile = ''
    names_listfile = ''
    seqids = []

    seqidre = regex.compile("\([a-zA-Z_0-9]+:[a-zA-Z0-9_\-]+\) *" +\
        "\([Begin:0-9]*\) *\([End:0-9]*\)")

    # Get list of seqIDs from input

    if string.find(input,'..') > -1:
        # Case of list file input
        lines = string.split(input,'\n')
        for line in lines:
            if seqidre.match(line) > -1:
                seqid = seqidre.group(1)
                rest = seqidre.group(2)
                seqids.append(seqid)
    else:
        # Case of seqID input
        seqids = [input]


    # Create a string that is the list file containing only seqIDs
    # (for GCG Names)
    names_listfile = "..\n" + string.joinfields(seqids,'\n')

    # Create file containing list file
    filename = writeFile(names_listfile, config)

    #############################################################
    # Step 2:                                                   #
    # Copy list file to GCG server.                             #
    #############################################################

    # Create temporary filename for output of GCG Names
    tempfile.tempdir = config.lookup('GCG_TEMP_DIR')
    cbrout = tempfile.mktemp()

    # Stamp time
    profiler.stamp('Before call to copyListFile')

    # Copy list file to GCG server
    copyListFile (filename,config.lookup('GCG_SERVER'),\
        config.lookup('GCG_ACCOUNT'),config.lookup('GCG_PASS'))

    #############################################################
    # Step 3:                                                   #
    # Construct input list file for GCG Names.                  #
    #############################################################

    # Stamp time
    profiler.stamp('Before call to runNames')

    # Run GCG Names and direct output to cbrout
    runNames (filename,cbrout,config.lookup('GCG_SERVER'),\
        config.lookup('GCG_ACCOUNT'))

    # Stamp time
    profiler.stamp('After call to runNames')

    # Return filename of GCG Names output file on remote GCG server
    return cbrout,profiler

def generateGCGListFile(tablist):

    # Purpose: Convert tab-delimited text file containing sequences
    #          into a GCG formatted list file.
    # Returns: A GCG formatted list file.
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    gcglist = "..\n"

    tablistre = regex.compile("\([A-Za-z \-]+\)\t" +\
        "\([A-Za-z0-9_\.\-]+\)\t*\([0-9]*\)\t*\([0-9]*\)")

    for line in string.split(tablist,'\n'):
        if tablistre.match(line) > -1:

            gcgline = ""
            begin = 0
            coorend = 0

            id_db = tablistre.group(1)
            id = tablistre.group(2)

            # Convert logicalDB names (names of sequence providers)
            # and seqIDs to GCG sequence specification (GCG database 
            # logical with seqID).

            if id_db == "Sequence DB":
                gcgline = "gball:%s" % id
            elif id_db == "RefSeq":
                gcgline = "refseqall:%s" % id
            elif id_db == "SWISS-PROT" or id_db == "TrEMBL":
                gcgline = "sptrnrdb:%s" % id
            elif id_db == "TIGR Mouse Gene Index":
                gcgline = "tigrmgi:%s" % id
            elif id_db == "TIGR Rat Gene Index":
                gcgline = "tigrrgi:%s" % id
            elif id_db == "TIGR Human Gene Index":
                gcgline = "tigrhgi:%s" % id
            elif id_db == "DoTS Mouse":
                # Translate DoTS IDs from 'dt.55133146' to 'dt55133146'
                id = regsub.gsub('\.','',id)
                gcgline = "dotsm:%s" % id
            elif id_db == "NIA Mouse Gene Index":
                gcgline = "niamgi:%s" % id

            # Add coordinates of subsequences if specified in input.

            if tablistre.group(3):
                begin = tablistre.group(3)
            if tablistre.group(4):
                coorend = tablistre.group(4)

            if begin > 0 and coorend > 0:
                gcgline = gcgline + "   Begin:%s   End:%s" %\
                    (begin,coorend)
            elif begin > 0:
                gcgline = gcgline + "   Begin:%s   End:" % (begin)
            elif coorend > 0:
                gcgline = gcgline + "   Begin:   End:%s" % (coorend)

            gcglist = gcglist + gcgline + "\n"

    return gcglist

#####################
# Private Functions #
#####################

def printString(s):

    # Purpose: Used as the callback function when using ftplib.retrlines()
    #          function.  It keeps appending the input string, s, to the
    #          global variable, contents.
    # Returns: nothing
    # Assumes: The global variable contents has been initialized.
    # Effects: nothing
    # Throws:  nothing

    global contents

    contents = contents + s + "\n"

def runNames (listfile,filename,gcgserver,login):

    # Purpose: Run the GCG Names program.
    # Returns: A GCG formatted list file.
    # Assumes: The input file list file resides on the remote
    #          GCG server.
    # Effects: nothing
    # Throws:  nothing

    # Construct command (using rsh)
    cmd = 'rsh -n -l %s %s names -in=@%s -def -out=%s' %\
        (login,gcgserver,listfile,filename)

    stdoutput, stderror, exitcode = runCommand.runCommand (cmd, {})


def runToFASTA (cbrout, filename,gcgserver,login):

    # Purpose: Run the GCG ToFASTA program.
    # Returns: A file containing FASTA-formatted sequences.
    # Assumes: The input file list file resides on the remote
    #          GCG server.
    # Effects: nothing
    # Throws:  nothing

    # Construct command (using rsh)
    cmd = 'rsh -n -l %s %s tofasta @%s -def -out=%s' %\
        (login,gcgserver,filename,cbrout)

    stdoutput, stderror, exitcode = runCommand.runCommand (cmd, {})


def getGCGOutput (filename,gcgserver,login,password):

    # Purpose: Get contents of a file from remote GCG server using FTP.
    # Returns: A string containing the file contents.
    # Assumes: nothing
    # Effects: Removes file that is read.
    # Throws:  error if unable to get file contents.

    global contents

    try:

        # Use FTP to transfer output of ToFASTA from remote GCG server
        contents = ''

        ftp = ftplib.FTP(gcgserver)
        ftp.login(login,password)
        cmd = 'RETR %s' % (filename)
        ftp.retrlines(cmd,printString)
        ftp.quit()

        stdoutput = contents

        # Need to add code back in that would remove filename
        # file on remote GCG server.

    except:
        raise error, 'Cannot get file from GCG server: %s' % filename

    return stdoutput
    

def copyListFile (filename, gcgserver, login, password):

    # Purpose: Put a file on a remote server.
    # Returns: nothing
    # Assumes: nothing
    # Effects: Creates a new file on a remote server.
    # Throws:  error if unable to put file on remote server.

    # Use FTP to transfer list of sequences to remote GCG server
    try:
        ftp = ftplib.FTP(gcgserver)
        ftp.login(login,password)
        cmd = 'STOR %s' % (filename)
        ftp.storlines(cmd,open (filename,'r'))
        ftp.quit()
    except:
        raise error, 'Cannot transfer file %s to GCG server' % filename

def seqReader(loopseq,config):

    # Purpose: Get the complete sequence for a given seqID from a
    #          string containing some number of FASTA-formatted
    #          sequences regardless of whether the sequence has
    #          been split by GCG.  Thus, all split sequences are
    #          joined by this function.
    # Returns: Three strings:
    #            1. base seqID of sequence that will be returned
    #            2. base seqID of next sequence in loopseq
    #         3. FASTA-format sequence
    #         4. remainder of loopseq
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    #############################################################
    # Algorithm:                                                #
    #                                                           #
    # Read first line to determine whether the sequence has     #
    # been split or not.  If not, then read through input until #
    # a new description line is found.  If split, then read     #
    # through all parts of the split sequence until it finds a  #
    # description line with a different base seqID or a the     #
    # first part of a split sequence (e.g., AC134402_0).        #
    #                                                           #
    # The different scenarios of sequences that may be          #
    # encountered in loopseq are:                               #
    # 1. nonsplitseqID to different nonsplitseqID with          #
    #    different baseseqID                                    #
    #    1. save buffer to seq, flag = 0                        #
    #    2. return seq                                          #
    # 2. nonsplitseqID to same nonsplitseqID with same          #
    #    baseseqID (duplicate sequences)                        #
    #    1. save buffer to seq, flag = 0                        #
    #    2. return seq                                          #
    # 3. nonsplitseqID to first splitseqID with different       #
    #    baseseqID                                              #
    #    1. save buffer to seq, flag = 0                        #
    #    2. return seq                                          #
    # 4. any splitseqID to next splitseqID with same baseseqID  #
    #    1. append buffer to combinebuffer                      #
    # 5. last splitseqID to first splitseqID with different     #
    #    baseseqID                                              #
    #    1. append buffer to combinebuffer                      #
    #    2. save combinebuffer to seq, flag = 1                 #
    #    3. return seq                                          #
    # 6. last splitseqID to first splitseqID with same          #
    #    baseseqID (duplicate sequences)                        #
    #    1. append buffer to combinebuffer                      #
    #    2. save combinebuffer to seq, flag = 1                 #
    #    3. return seq                                          #
    # 7. last splitseqID to nonsplitseqID with different        #
    #    baseseqID                                              #
    #    1. append buffer to combinebuffer                      #
    #    2. save combinebuffer to seq, flag = 1                 #
    #    3. return seq                                          #
    # 8. nonsplitseqID to END OF FILE                           #
    #    1. save buffer to seq, flag = 0                        #
    #    2. return seq                                          #
    # 9. last splitseqID to END OF FILE                         #
    #    1. append buffer to combinebuffer                      #
    #    2. save combinebuffer to seq, flag = 1                 #
    #    3. return seq                                          #
    #############################################################

    buffer = ''
    combinebuffer = ''
    baseseqID = ''
    firstbaseseqID = ''
    nextbaseseqID = ''
    flag = 0
    newloopseq = ''
    lastindex = 0
    firstsplitflag = 0
    nextsplitflag = 0

    # Regular expression to find first line of a sequence in a 
    # multi-sequence FASTA-format file, and capture seqID.
    seqidre = regex.compile(">\([a-zA-Z0-9_\-]+\) .*")

    # Regular expression to determine whether a sequence 
    # has been split by GCG, and get the base seqID.  For example, 
    # the base seqID of 'AC134402_0' is 'AC134402'.
    splitseqidre = regex.compile("\([a-zA-Z0-9][a-zA-Z0-9\-]_*" +\
        "[a-zA-Z0-9\-]+\)_\([0-9]+\)")

    # Get first sequence (seqID and base seqID from first line)
    if seqidre.match(loopseq) > -1:
        firstseqID = seqidre.group(1)
        if splitseqidre.match(firstseqID) > -1:
            firstsplitflag = 1
            firstbaseseqID = splitseqidre.group(1)
            combinebuffer = ''
        else:
            firstbaseseqID = firstseqID

    # Save the firstbaseseqID so that it can be returned
    savefirstbaseseqID = firstbaseseqID

    # Parse all sequences except for last one
    while (lastindex > -1):
        # Find next description line
        lastindex = string.find(loopseq,'>',string.find(loopseq,\
            '\n')+1)

        if lastindex > -1:
            # Store previous sequence in buffer
            buffer = string.rstrip(loopseq[:lastindex])

            # Reset loopseq so it doesn't include previous sequence
            loopseq = loopseq[lastindex:]
        
            # If description line then get baseseqID of next 
            # sequence
            if seqidre.match(loopseq) > -1:
                nextseqID = seqidre.group(1)
                if splitseqidre.match(nextseqID) > -1:
                    nextbaseseqID = splitseqidre.group(1)
                    splitid = splitseqidre.group(2)
                    nextsplitflag = 1
                else:
                    nextbaseseqID = nextseqID
                    nextsplitflag = 0

            # If new sequence is not a split sequence, then
            # handle scenarios #1, #2, #7
            if nextsplitflag == 0:
                # handle scenarios #1, #2
                if firstsplitflag == 0:
                    seq = buffer
                    #firstbaseseqID = nextbaseseqID
                    flag = 0
                    break

                # handle scenario #7
                else:
                    combinebuffer = \
                        fastaseqlib.appendSequences(combinebuffer,\
                        buffer, config.lookup('OVERLAP_LENGTH'),\
                        config.lookup('LINE_LENGTH'))
                    seq = combinebuffer
                    firstbaseseqID = nextbaseseqID
                    flag = 1
                    break
            else:
                # handle scenarios #3 and #5
                if nextbaseseqID != firstbaseseqID:
                    # handle scenario #3
                    if splitid != '0':
                        seq = buffer
                        flag = 0
                        firstbaseseqID = nextbaseseqID
                        break

                    # handle scenario #5
                    else:
                        combinebuffer = \
                            fastaseqlib.appendSequences(combinebuffer,\
                            buffer, config.lookup('OVERLAP_LENGTH'),\
                            config.lookup('LINE_LENGTH'))
                        seq = combinebuffer
                        firstbaseseqID = nextbaseseqID
                        flag = 1
                        break
        
                # handle scenario #4 and #6
                else:
                    # handle scenario #4
                    if splitid != '0':
                        combinebuffer = \
                            fastaseqlib.appendSequences(combinebuffer,\
                            buffer, config.lookup('OVERLAP_LENGTH'),\
                            config.lookup('LINE_LENGTH'))
                        firstbaseseqID = nextbaseseqID

                    # handle scenario #6
                    else:
                        combinebuffer = \
                            fastaseqlib.appendSequences(combinebuffer,\
                            buffer, config.lookup('OVERLAP_LENGTH'),\
                            config.lookup('LINE_LENGTH'))
                        seq = combinebuffer
                        firstbaseseqID = nextbaseseqID
                        flag = 1
                        break

        # Parse last sequence
        else:

            # Store previous sequence in buffer
            buffer = loopseq

            loopseq = ''

            # handle scenario #8
            if firstsplitflag == 0:
                seq = buffer
                flag = 0

            # handle scenario #9
            else:
                combinebuffer = \
                    fastaseqlib.appendSequences(combinebuffer,\
                    buffer, config.lookup('OVERLAP_LENGTH'),\
                    config.lookup('LINE_LENGTH'))
                seq = combinebuffer
                flag = 1        

            break

    newloopseq = loopseq

    return savefirstbaseseqID,nextbaseseqID,seq,newloopseq    


def mergeListFiles(nameslistfile,listfile):

    # Purpose: Merge two GCG formatted list files:
    #            1. (listfile) One that may contain subsequence
    #               coordinates, but no split sequences
    #               or no new versions of sequences.
    #            2. (nameslistfile) Another that is from the output 
    #               of GCG Names when run with the above list file
    #               as input.  This list file does not include 
    #               coordinates, but may contain split sequences 
    #               and/or new versions of sequences.
    #          Merge the above list files into one where:
    #            1. any sequence in the first file that is split 
    #               up by GCG, is replaced with the identifiers
    #               of all the split subsequences.
    #            2. any sequence in the first file that has a 
    #               newer version than in the release is replaced
    #               with that newer version only.
    # Returns: 1. string containing merged GCG format list file
    #          2. list of seqids in first list file
    #          3. dictionary of split subsequence coordinates
    #             keyed by baseseqID
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    #############################################################
    # Algorithm:                                                #
    #                                                           #
    # 1. Parse first list file to get:                          #
    #      1. logical DB IDs                                    #
    #      2. seqIDs                                            #
    #      3. subsequence coordinates (if any)                  #
    # 2. Parse second list file to get:                         #
    #      1. logical DB IDs                                    #
    #      2. seqIDs                                            #
    #      3. base seqIDs for any split sequences               #
    # 3. Determine which sequences in second list file:         #
    #      1. are split up into multiple sequences              #
    #      2. have a new version                                #
    # 4. Construct a list of logical DB IDs and seqIDs          #
    #    for split and non-split sequences such that only       #
    #    new versions of the sequences are included.            #
    # 5. Construct GCG format list file using above list and    #
    #    add back any subsequence coordinates specified in      #
    #    the original list file for only non-split sequences.   #
    #    For split sequences, any subsequence coordinates       #
    #    are saved in a dictionary that is returned by the      #
    #    function.                                              #
    #############################################################

    splitsubseqs = {}
    seqids = []
    coords = {}
    mergedlistfile = "..\n"
    gcgnamesseqids = []
    gb_news = {}
    namesseqs = {}
    baseseqids = {}
    newgcgnamesseqids = []
    missingseqs = []

    # Regex to get logical DB ID and seqID as one string, along with 
    # begin and end coordinates as separate strings.
    linere = regex.compile("\([A-Za-z_0-9]+:[a-zA-Z0-9_\-]+\) *" +\
        "[Begin:]*\([0-9]*\) *[End:]*\([0-9]*\)")

    # Regex to get logical DB ID and seqID as separate strings, and
    # begin and end coordinates as separate strings.
    lineseqidre = regex.compile("\([A-Za-z_0-9]+\):\([a-zA-Z0-9_\-]+\) *" +\
        "[Begin:]*\([0-9]*\) *[End:]*\([0-9]*\)")

    # Regex to get logical DB ID and seqID as separate strings
    seqidre = regex.compile("\([A-Za-z_0-9]+\):\([a-zA-Z0-9_\-]+\)")

    # Regex to detect split sequence ID and get base seqID
    splitseqidre = regex.compile("\([a-zA-Z0-9][a-zA-Z0-9\-]_*" +\
        "[a-zA-Z0-9]+\)_[0-9]+")

    # Regex to get coordinates from a list
    coordre = regex.compile("\([0-9]*\),\([0-9]*\)")

    # Process contents of original list file to get coordinates

    listfilelines = string.split(listfile,'\n')

    for line in listfilelines:

        begin = 0
        coorend = 0

        # Determine whether line is correct format
        if lineseqidre.match(line) > -1:

            seqid = lineseqidre.group(2)
            # Translate lowercase to uppercase
            seqid = string.upper(seqid)

            seqids.append(seqid)

            # Get begin and end subsequence coordinates if any
            if lineseqidre.group(3):
                begin = lineseqidre.group(3)

            if lineseqidre.group(4):
                coorend = lineseqidre.group(4)
            
            # Store subsequence coordinates

            if begin > 0 and coorend > 0:
                if coords.has_key(seqid):
                    coords[seqid] = coords[seqid] +\
                        ["%s,%s" % (begin,coorend)]
                else:
                    coords[seqid] = ["%s,%s" % (begin,coorend)]
            elif begin > 0:
                if coords.has_key(seqid):
                    coords[seqid] = coords[seqid] + ["%s," % (begin)]
                else:
                    coords[seqid] = ["%s," % (begin)]
            elif coorend > 0:
                if coords.has_key(seqid):
                    coords[seqid] = coords[seqid] + [",%s" % (coorend)]
                else:
                    coords[seqid] = [",%s" % (coorend)]


    # Process second list file

    namelistlines = string.split(nameslistfile,'\n')
    for line in namelistlines:
        if linere.match(line) > -1:
            gcgnamesseqid = linere.group(1)
            # Translate lowercase to uppercase
            gcgnamesseqid = string.upper(gcgnamesseqid)
            gcgnamesseqids.append(gcgnamesseqid)

            if seqidre.match(gcgnamesseqid) > -1:
                logicaldb = seqidre.group(1)
                dbid = seqidre.group(2)
                if splitseqidre.match(dbid) > -1:

                    # If a split sequence, store baseseqid

                    baseseqid = splitseqidre.group(1)
                    if not baseseqids.has_key(baseseqid):
                        baseseqids[baseseqid] = baseseqid
                    
                    # Store info if sequence is split and from 'gb_new' 

                    if logicaldb == 'GB_NEW':
                        if not gb_news.has_key(baseseqid):
                            gb_news[baseseqid] = {}
                            gb_news[baseseqid][logicaldb + ':' +\
                                dbid] = logicaldb + ':' + dbid
                        else:
                            gb_news[baseseqid][logicaldb + ':' +\
                                dbid] = logicaldb + ':' + dbid

                    # Store info for other split sequences

                    else:
                        if not namesseqs.has_key(baseseqid):
                            namesseqs[baseseqid] = {}
                            namesseqs[baseseqid][logicaldb + ':' +\
                                dbid] = logicaldb + ':' + dbid
                        else:
                            namesseqs[baseseqid][logicaldb + ':' +\
                                dbid] = logicaldb + ':' + dbid

                # Store info for non-split sequences

                else:
                    if logicaldb == 'GB_NEW':
                        if not gb_news.has_key(dbid):
                            gb_news[dbid] = {}
                            gb_news[dbid][logicaldb + ':' + dbid] =\
                                logicaldb + ':' + dbid
                    else:
                        if not namesseqs.has_key(dbid):
                            namesseqs[dbid] = {}
                            namesseqs[dbid][logicaldb + ':' + dbid] =\
                                logicaldb + ':' + dbid


    # Process list of seqIDs in second list file to make sure
    # only newest copy of a sequence is specified.

    for seqid in seqids:
        keys = []

        # Check to see whether seqid has been split in names list
        # file
        if baseseqids.has_key(seqid):

            # Check to see whether seqid in 'gb_new'
            if gb_news.has_key(seqid):

                # Construct merged list file record for sequence
                # even if it has been split up into chunks
                keys = gb_news[seqid].keys()
                keys.sort()

                for id in keys:
                    newgcgnamesseqids.append(id)
            else:

                if namesseqs.has_key(seqid):
                    # Construct merged list file record for sequence
                    # even if it has been split up into chunks
                    keys = namesseqs[seqid].keys()
                    keys.sort()

                    for id in keys:
                        newgcgnamesseqids.append(id)
                else:
#                    raise error, 'Sequence %s not found' % seqid
                    missingseqs.append(seqid)

        # Process sequences that have not been split
        else:

            # Check to see whether seqid in 'gb_new'
            if gb_news.has_key(seqid):

                # Construct merged list file record for sequence
                # even if it has been split up into chunks
                keys = gb_news[seqid].keys()
                keys.sort()

                for id in keys:
                    newgcgnamesseqids.append(id)
            else:

                if namesseqs.has_key(seqid):
                    # Construct merged list file record for sequence
                    # even if it has been split up into chunks
                    keys = namesseqs[seqid].keys()
                    keys.sort()

                    for id in keys:
                        newgcgnamesseqids.append(id)
                else:
                    missingseqs.append(seqid)

    gcgnamesseqids = newgcgnamesseqids

    # Generate merged list file

    for gcgnamesseqid in gcgnamesseqids:

        if seqidre.match(gcgnamesseqid) > -1:
            seqid = seqidre.group(2)

        # Process split sequences

        if splitseqidre.match(seqid) > -1:
            baseseqid = splitseqidre.group(1)

            # Save coordinates to splitsubseqs dictionary

            if coords.has_key(baseseqid):
                splitsubseqs[baseseqid] = coords[baseseqid]

            # Make entry in merged list file

            mergedlistfile = mergedlistfile + gcgnamesseqid + "\n"
    
        # Process unsplit sequences

        else:

            # If subsequence coordinates

            if coords.has_key(seqid):
                coord = coords[seqid][0]
                coords[seqid] = coords[seqid][1:]
                begin = 0
                coorend = 0

                # Append subsequence coordinates if needed.
                # Could also use string.split instead of regex.

                if coordre.match(coord) > -1:
                    if coordre.group(1):
                        begin = coordre.group(1)
                    if coordre.group(2):
                        coorend = coordre.group(2)

                    # Append line to merged list file appropriately

                    if begin > 0 and coorend > 0:
                        mergedlistfile = mergedlistfile + \
                            "%s   Begin:%s  End:%s\n" %\
                            (gcgnamesseqid,begin,coorend)
                    elif begin > 0:
                        mergedlistfile = mergedlistfile + \
                            "%s   Begin:%s  End:\n" %\
                            (gcgnamesseqid,begin)
                    elif coorend > 0:
                        mergedlistfile = mergedlistfile + \
                            "%s   Begin:  End:%s\n" %\
                            (gcgnamesseqid,coorend)

            # If no subsequence coordinates

            else:
                mergedlistfile = mergedlistfile + gcgnamesseqid + "\n"

    return mergedlistfile,seqids,splitsubseqs,\
        string.joinfields(missingseqs,',')+"\n"
    
def writeFile(content, config):

    # Purpose: Writes a string to a file on the local file system.
    # Returns: A string containing the name of the new file.
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    tempfile.tempdir = config.lookup('LOCAL_TEMP_DIR')
    filename = tempfile.mktemp()
    fd = open (filename, 'w')
    fd.write (content)
    fd.close()

    return filename


def processSequences (sequence, seqids, splitsubseqs, config):

    # Purpose: Processes a group of sequences to make sure that:
    #          1. All sequences from original list file in output.
    #          2. Any specified subsequences from sequences that
    #             are split up by GCG are returned.
    # Returns: A string containing a group of sequences in 
    #          FASTA format.
    # Assumes: nothing
    # Effects: nothing
    # Throws:  nothing

    #############################################################
    # ALGORITHM:                                                #
    #     1. Get the first sequence from the ToFASTA output     #
    #        and check to see that it is the first sequence     #
    #        in the original list file.  This sequence will     #
    #        be either be the complete sequence (regardless of  #
    #        whether it was split up by GCG), or a subsequence  #
    #        of a non-split up sequence if coordinates          #
    #        specified in original list file.  If the sequence  #
    #        was split up by GCG and the coordinates for that   #
    #        sequence was specified in the original list file,  #
    #        then the subsequence is generated. The first       #
    #        sequence (or subsequence) is saved in the variable #
    #        'store'.                                           #
    #     2. Get the next sequence in the remaining ToFASTA     #
    #        output like what was done with the previous        #
    #        sequence.  If this sequence is the next one in the #
    #        original list file, then the variable 'store' is   #
    #        appended to 'masterbuffer', and the variable       #
    #        'store' is reset with the current sequence.        #
    #        Note that this will allow for removing duplicate   #
    #        sequences returned from ToFASTA.                   #
    #     3. All remaining sequences are processed as described #
    #        above.                                             #
    #############################################################

    # Set flag to tell whether a seq has been split in GCG
    flag = 0

    # Dictionary of seqIDs found in ToFASTA output
    foundseqids = {}

    # String containing a buffer that is about to be saved.
    # 'store' is used to keep the most recent copy a given 
    # sequence that may have multiple copies in the ToFASTA
    # output.
    store = ''

    # String containing the output sequences as they are processed
    masterbuffer = ''

    # String containing the portion of ToFASTA output not yet
    # processed.
    loopsequence = sequence

    # Scenarios of different sequences:
    # 1. current and next baseseqIDs are different:
    #    1. save first sequence to masterbuffer
    #    2. save next sequence to store
    # 2. current and next baseseqIDs are same (duplicate sequences)
    #    1. save next sequence to store
    # 3. if final sequence
    #    1. save final sequence to masterbuffer


    # Read and store sequences    
    loop_counter = 0
    while (len(loopsequence) > 0):

        currentbaseseqID,nextbaseseqID,nextseq,loopsequence = \
            seqReader(loopsequence,config)

        # Check to see if sequence has been split by GCG and has a 
        # subsequence specified in the original list file.
        if splitsubseqs.has_key(currentbaseseqID):
            subseqcoords = string.split(\
                splitsubseqs[currentbaseseqID][0],',')
            splitsubseqs[currentbaseseqID] = \
                splitsubseqs[currentbaseseqID][1:]
            nextseq = fastaseqlib.getSubSeq(nextseq,subseqcoords[0],\
                subseqcoords[1],config)
            if loop_counter == 0:
                store = nextseq = "\n"

        # handle scenario #1
        elif (currentbaseseqID != nextbaseseqID) and (currentbaseseqID == seqids[0]):
            if loop_counter == 0:
                store = nextseq + "\n"
            else:
                masterbuffer = masterbuffer + store
                store = nextseq + "\n"

            # Adjust seqids to remove missing ID(s)
            id_loop_counter = 0
            for id_loop in seqids:
                if currentbaseseqID == id_loop:
                    seqids = seqids[id_loop_counter+1:]
                    break
                id_loop_counter = id_loop_counter + 1

        # handle scenario where current sequence does not match next
        # sequence in original list file
        elif (currentbaseseqID != seqids[0]):
            if loop_counter == 0:
                store = nextseq + "\n"
            else:
                masterbuffer = masterbuffer + store
                store = nextseq + "\n"
            # Adjust seqids to remove missing ID(s)
            id_loop_counter = 0
            for id_loop in seqids:
                if currentbaseseqID == id_loop:
                    seqids = seqids[id_loop_counter+1:]
                    break
                id_loop_counter = id_loop_counter + 1


        loop_counter = loop_counter + 1

    # handle scenario #3
    masterbuffer = masterbuffer + store

    sequence = masterbuffer

    return sequence



