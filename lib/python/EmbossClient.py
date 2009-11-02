#!./python

# Module: EmbossClient.py
# Author: marka
#
# Purpose: This object is designed to act as a web service client to the 
#	soaplab2 webservice wrapper for EMBOSS.  EmbossClient objects must be
#	instantiated with a seqfetch config object.  The only public method is 
#	getSequences which takes a string containing newline separated  query 
#	strings received by tofasta.cgi.  This string is parsed into formatted 
#	query strings which are inserted into a webservice soap request to soaplab.
#	Each requested sequence can generate multiple requests.  This is because 
#	each requested sequence id may require queries to multiple EMBOSS 
#	databases.  This is currently true for GenBank and RefSeq sequence ids.  
#	The response from soaplab is then parsed for a fault and the sequence.
#	A string containing the requested sequences in fasta format AND a string  
#	containing a list of the sequences not found are returned by getSequences.
#	Some possible error conditions that will raise an error include; httplib
#	errors connecting to soaplab web service, return of a fault by web service,
#	or no sequences parsed from input list.
#
# Notes: 
#	1. This object requires a seqfetch config object at instantiation.  This 
#	is mainly due to the fact that the logginf methods provided by 
#	tofastalib.py require a config object to obtain the logfile path.  
#	2. This object requires lxml for xml processing.
#	3. Example request namespaces removed for brevity:
#			<SOAP-ENV:Envelope xmlns:org="http:..." 
#					xmlns:encodingStyle="http:..." 
#					xmlns:SOAP-ENV="http:...>
#				<SOAP-ENV:Header/>
#				<SOAP-ENV:Body>
#					<org:runAndWaitFor>
#						<arg0>
#							<entries>
#								<key>sequence_usa</key>
#								<value>
#									<singleString>gbnew:AY174014</singleString>
#								</value>
#							</entries>
#						</arg0>
#					</org:runAndWaitFor>
#				</SOAP-ENV:Body>
#			</SOAP-ENV:Envelope>
#	4. Example response namespaces removed and element values for brevity:
#		<Envelope>
#			<Body>
#				<runAndWaitForResponse>
#				<return>
#					<entries>
#						<key>outseq_url</key>
#						<value><singleString>http:...</singleString></value>
#					</entries>
#					<entries>
#						<key>detailed_status</key>
#						<value><singleString>0</singleString></value>
#					</entries>
#					<entries>
#						<key>report</key>
#						<value><singleString>Summary:.../singleString></value>
#					</entries>
#					<entries>
#						<key>outseq</key>
#						<value><singleString>&gt;X...</singleString></value>
#					</entries>
#				</return>
#				<runAndWaitForResponse>
#			<Body>
#		<Envelope>


import sys, httplib, re, string
from types import *
import lxml.etree
import tofastalib

class EmbossClient:
	"""
	This class acts as a web services client to the soaplab2 installation
	on hygelac that exposes the EMBOSS commands as web services.  
	
	A valid config object is required on order for this class to get the 
	necessary EMBOSS_HOST and EMBOSS_URL parameter values.  Once 
	instantiated getSequences can be called with a list of sequence ids to
	fetch from EMBOSS
	"""
	
	# error string template for missing sequences message
	missing_seq = "Sequence %s not found"
	
	# error identifier
	error = "EmbossClient.py error"
	
	# dictionary to define query strings for each sequence type.
	# the key is the sequence type, and the value is a list of query 
	# strings that the id will be inserted into.  the query strings are    
	# in database:id format. queries will be run in the order defined.
	databases = {"Sequence DB":["gbnew:%s", "genbank:%s"], \
		"RefSeq":["refseqAll:%s"], \
		"SWISS-PROT":["sprot:%s"], \
		"TrEMBL":["tr:%s"], \
		"DFCI Mouse Gene Index":["mgi:%s"], \
		"DFCI Rat Gene Index":["rgi:%s"], \
		"DFCI Human Gene Index":["hgi:%s"], \
		"DoTS Mouse":["dots:%s"], \
		"NIA Mouse Gene Index":["nia:%s"]}
	
	def __init__(self, config):
		"""
		Constructor for EmbossClient requires a seqfetch config object
		"""
		# Purpose: returns 's' if it is non-blank, or a defined
		#	string (self.blank) to substitute for a blank 's'
		# Returns: string
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing
		
		if config:
			self.config = config
		
		return
	
	
	def __buildRequest(self, 
				seq	# a list of query strings to insert into requests
				):
		"""
		Purpose:	Builds a SOAP request for the given seq string
		Returns:	A list of strings containing the soap requests for the
					given input query strings
		Assumes:	The seq contains database:id formatted strings
		Effects:	Nothing
		Throws:		Nothing
		"""
		
		# define xml namespaces used		
		ec = "http://schemas.xmlsoap.org/soap/encoding/"
		soapEnv = "http://schemas.xmlsoap.org/soap/envelope/"
		org="http://org.soaplab.2"
		
		# setup namespace strings
		SOAP = "{%s}" % soapEnv
		ORG = "{%s}" % org
		EC = "{%s}" % ec
		
		# create namespace dictionary used to build soap request
		NSMAP = {'SOAP-ENV':soapEnv, 'org':org, 'encodingStyle':ec}
		
		# list to hold generated soap requests 
		requests = []
		
		# build soap request for each sequence query string.  
		# see example request in module comments.
		for s in seq:
			# setup soap request elements: envelope, header, body
			domdoc = lxml.etree.Element(SOAP + "Envelope", nsmap=NSMAP)
			header = lxml.etree.SubElement(domdoc, SOAP + "Header")
			body = lxml.etree.SubElement(domdoc, SOAP + "Body")
			# setup soaplab request elements
			run = lxml.etree.SubElement(body, ORG + "runAndWaitFor")
			arg0 = lxml.etree.SubElement(run, "arg0")
			entries = lxml.etree.SubElement(arg0, "entries")
			key = lxml.etree.SubElement(entries, "key")
			# define sequence format as Uniform Sequence Address
			key.text = "sequence_usa"
			value = lxml.etree.SubElement(entries, "value")
			singleString = lxml.etree.SubElement(value, "singleString")
			# insert sequence query string
			singleString.text = s
			
			# debug log soap request
			tofastalib.writeToErrorLogDebug (lxml.etree.tostring(domdoc), self.config)
			
			# add request to list
			requests.append(lxml.etree.tostring(domdoc))
		
		# return request list
		return requests
	
	
	def __parseSequenceInput(self, 
				tablist  # a string containing newline separated seq requests
				):
		"""
		Purpose:	Parses the given list of sequence ids and translates them 
					a query string to be used in a SOAP request to soaplab
		Returns:	A list of strings containing the query strings for the 
					given list of sequence ids
		Assumes:	tablist contains a newline separated list of the query 
					strings received by tofasta.cgi
		Effects:	Nothing
		Throws:		Nothing
		"""
		
		list = [] # list to hold sequence query strings parsed from input
		
		# regular expression to capture sequence type, id, 
		# begin coordinates, and end coordinates
		tablistre = re.compile("([A-Za-z\s\-]+)\t" +\
			"([A-Za-z0-9_\.\-]+)\t*([0-9]*)\t*([0-9]*)")
		
		# parse list of request ids
		for line in string.split(tablist,'\n'):
			# get matches to regular expression
			match = tablistre.match(line)
			
			if match:				
				id_db = match.group(1)  # sequence type
				id = match.group(2) 	# sequence id
				
				# lookup sequence type to get corresponding query strings
				# insert id into query string and store in a list
				if id_db in self.databases:
					seq = []   # list to hold formatted query strings					
					for db in self.databases[id_db]:						
						seq.append(db % id)	
					list.append(seq)
				
				else:
					# no query string definition for sequence type, log debug
					msg = 'No query string definition for sequence type: %s'
					tofastalib.writeToErrorLogDebug (msg % id_db, self.config)
					
			elif line != '' and line != '..':
					# no query string definition for sequence type, log error
					# this may be normal, but good to check what falls through
					msg = 'No sequence data found in request string: %s'
					tofastalib.writeToErrorLog (msg % line, self.config)
		
		# raise error if no sequence ids detected
		if len(list) == 0:
			raise self.error, 'All sequences missing'
		
		# return list of formatted query strings
		return list
	
	
	def __embossClient(self, 
				requests	# a list of strings containing query strings to /
					   	   	# issue to the soaplab web service
				):
		"""
		Purpose:	Issues web service client requests to the soaplab seqret
					web service and returns the sequence or None
		Returns:	A string containing the requested sequence, or None
		Assumes:	EmbossClient contains a valid seqfetch config object that
					contains the EMBOSS_HOST and EMBOSS_URL parameters 
		Effects:	Nothing
		Throws:		Error if any httplib.HTTP error occurs
		"""
		
		try:
			# perform a web service queries to emboss for each request until
			# a sequence is returned
			for req in requests:
				# construct the header and send the request
				webservice = httplib.HTTP( self.config.lookup('EMBOSS_HOST') )
				webservice.putrequest( "POST", self.config.lookup('EMBOSS_URL') )
				webservice.putheader( "Host", self.config.lookup('EMBOSS_HOST') )
				webservice.putheader("Content-type","text/xml; charset=\"UTF-8\"")
				webservice.putheader( "Content-length", "%d" % len(req) )
				webservice.putheader( "SOAPAction", "org:runAndWaitFor" )
				webservice.endheaders()
				webservice.send(req)
				
				# get the response from the web service and parse the reply
				statuscode, statusmessage, header = webservice.getreply()
				f = webservice.getfile()
				seq = self.__parseReply(lxml.etree.parse(f))
				f.close()
				
				# return first sequence found
				if seq:
					return seq
		# communication error, log error and raise to caller
		except Exception, inst:
			tofastalib.writeToErrorLog (inst, self.config)
			raise self.error, 'Unable to connect to EMBOSS server.'
		
		return
	
	
	def __parseReply(self, 
					reply  # 
					):
		"""
		Purpose:	Parses a web service reply.  Evaluates reply for a fault,
					parses sequence string from reply, logs error conditions.
		Returns:	A string containing the contents of the outseq element in
					reply, or None if outseq is absent from reply.
		Assumes:	Replay is an xml element object
		Effects:	Nothing
		Throws:		Error if a web service fault is in reply
		"""
		
		# get faultstring element
		fault = reply.xpath('//faultstring')
		# get <singleString> element that is a child of <entries> and has a <key> 
		# sibling with the value 'report'
		report = reply.xpath('//entries[key="report"]/value/singleString')
		# get <singleString> element that is a child of <entries> and has a <key> 
		# sibling with the value 'outseq'
		outseq = reply.xpath('//entries[key="outseq"]/value/singleString')

		if len(fault) < 1:
			# return sequence(s) from reply
			if len(outseq) > 0:
				for seq in outseq:
					if seq is not None:
						return seq.text
			else:
				# no sequence found, return None
				
				# enable if you want to see soaplab output
				# tofastalib.writeToErrorLogDebug (report[0].text, self.config)
				return None
		else:
			# soaplab returned a fault, log and raise error
			message = 'EMBOSS Server error: ' + fault[0].text
		  	tofastalib.writeToErrorLogDebug (message, self.config)
		  	raise self.error, message
	
	
	def getSequences(self, 
					seqs   # a string containing newline separated seq requests
					):
		"""
		Purpose:	Retrieve fasta sequence strings for the given list of 
					sequence ids.
		Returns:	A string containing the fasta sequences for the given 
					sequence ids and a string containing the list of ids not
					found in EMBOSS
		Assumes:	Reply is an xml element object
		Effects:	Nothing
		Throws:		Any error generated internally by EmbossClient
		"""
		
		replies = []  # list if sequences returned
		missing = []  # list of missing sequences
		
		try:
			# parse ids from input and retrieve sequences
			for request in self.__parseSequenceInput(seqs):
				if request:
					# build request and call embossClient to retrieve sequence
					result = self.__embossClient( self.__buildRequest(request) )
					# add result to appropriate list
					if result:
						replies.append( result )
					else:
						# sequence not found, log debug, add to missing list
						message = request[-1].split(':')[-1]
						tofastalib.writeToErrorLogDebug ('Sequence %s not found' % message, self.config)
						missing.append( message )
			
		except self.error, message:
			raise self.error, '' + message
		
		# return results
		return '\n%s' % '\n'.join(replies), '%s\n' % ','.join(missing)
