# Name:         CGI.py
# Purpose:      defines the CGI class (see below for details)

import time
import sys
import os
import cgi
import types

import errorlib
import mgi_cgi

class CGI:
        # Concept:
        #       IS:   a CGI script, complete with logging and exception
        #             handling
        #       HAS:  an error log and some parameters passed in
        #       DOES: parses parameters, logs messages, wraps main "program"
        #             (method) in error handling
        # Implementation:

        def __init__ (self,
                logfilename = None      # string; name of log file to create
                ):
                # Purpose: initializes the object and creates a log file with
                #       the specified name, if possible.  If no log file name
                #       is specified, we do not do logging.
                # Returns: nothing
                # Assumes: nothing
                # Effects: creates a log file in the OS
                # Throws: nothing
                # Notes: If we cannot create the log file, we write a message
                #       to sys.stderr.  So, look there if you are not finding
                #       the specified file.

                if logfilename is None:
                        self.logfd = None
                else:
                        try:
                                self.logfd = open (logfilename, 'w')
                        except:
                                self.logfd = sys.stderr
                                self.log ('Could not open log file %s' % \
                                        logfilename)
                                self.logfd = None
                self.fields = {}
                #print 'Content-type: text/html'
                #print
                return

        def log (self,
                message         # string; message to write to log file
                ):
                # Purpose: writes a message to the log file (if one exists)
                # Returns: nothing
                # Assumes: nothing
                # Effects: writes to the log file
                # Throws: nothing
                # Notes: The format of the line written includes:
                #       script name: [date & time] [client IP] message\n

                if self.logfd is not None:
                        cginame = os.path.basename (os.environ['SCRIPT_NAME'])
                        datetime = time.asctime (time.localtime (time.time()))
                        ip = os.environ['REMOTE_ADDR']
                        self.logfd.write ('%s: [%s] [client %s] %s\n' % \
                                (cginame, datetime, ip, message))
                return

        def get_parms (self,
                default_fields = None,  # dict; { fieldname : {
                                        #           'op' : default operator,
                                        #           'val' : default value } }
                default_types = None    # dict; { fieldname : type string }
                ):
                # Purpose: get parameters passed into CGI script by GET or
                #       POST method
                # Returns: see Notes
                # Assumes: nothing
                # Effects: reads from stdin
                # Throws: propagates exceptions from mgi_cgi.FieldStorage
                # Notes: The default behavior of this method is to return a
                #       dictionary which has fieldnames as keys.  Each field-
                #       name is mapped to either a string (for single-valued
                #       fields) or a list of strings (for multi-valued
                #       fields).
                #       If the user specifies default_fields and default_types
                #       we instead act as a wrapper for the mgi_cgi library's
                #       FieldStorage class, an instance of which we return.
                #       For examples of this behavior, see the WI, which uses
                #       the mgi_cgi library extensively.

                if not ((default_fields is None) and (default_types is None)):
                        self.fields = mgi_cgi.FieldStorage ( \
                                default_fields, default_types)
                else:
                        fs = cgi.FieldStorage()
                        for key in list(fs.keys()):
                                if type(fs[key]) == list:
                                        self.fields[key] = []
                                        for item in fs[key]:
                                                self.fields[key].append ( \
                                                        item.value)
                                else:
                                        self.fields[key] = fs[key].value
                return self.fields

        def go (self,
                handler = errorlib.handle_error # function to call if an
                                                # exception occurs
                ):
                # Purpose: wraps the main() method in exception handling
                # Returns: nothing
                # Assumes: nothing
                # Effects: runs the main() method and outputs an explanatory
                #       HTML page if an error occurs
                # Throws: nothing

                try:
                        self.main()
                except SystemExit:
                        pass
                except:
                        handler()


        def main (self):
                # Purpose: abstract method.  Conceptually, this is the "main
                #       program" of the CGI script.  Define this in a subclass
                #       for each CGI script.
                # Returns: nothing
                # Assumes: nothing
                # Effects: nothing
                # Throws: nothing

                return

#
# Warranty Disclaimer and Copyright Notice
# 
#  THE JACKSON LABORATORY MAKES NO REPRESENTATION ABOUT THE SUITABILITY OR 
#  ACCURACY OF THIS SOFTWARE OR DATA FOR ANY PURPOSE, AND MAKES NO WARRANTIES, 
#  EITHER EXPRESS OR IMPLIED, INCLUDING MERCHANTABILITY AND FITNESS FOR A 
#  PARTICULAR PURPOSE OR THAT THE USE OF THIS SOFTWARE OR DATA WILL NOT 
#  INFRINGE ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS, OR OTHER RIGHTS.  
#  THE SOFTWARE AND DATA ARE PROVIDED "AS IS".
# 
#  This software and data are provided to enhance knowledge and encourage 
#  progress in the scientific community and are to be used only for research 
#  and educational purposes.  Any reproduction or use for commercial purpose 
#  is prohibited without the prior express written permission of the Jackson 
#  Laboratory.
# 
# Copyright 1996, 1999, 2002 by The Jackson Laboratory
# All Rights Reserved
#
