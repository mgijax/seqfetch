#!/usr/local/bin/python

# Name: Profiler.py
# Purpose: to provide an easy and fast mechanism for software engineers to
#	do profiling of their code.  (Python's deterministic profiler,
#	while providing very complete results, is horribly slow.)  Use of
#	this Profiler class will allow the programmer to determine the level
#	of granularity needed, by adding "stamp()" calls through his/her
#	code as needed.
# Notes:  If run as a script, this module runs some code to give example
#	output from the Profiler.

import copy
import regsub
import sys
import time
import string

class Profiler:
	# IS: A Profiler is an object to help in measuring the efficiency of
	#	pieces of code.
	# HAS: an ordered list of time-stamped messages.  Time-stamping starts
	#	at 0.0 seconds when the Profiler is instantiated.
	# DOES: time-stamps and stores messages, writes those messages to
	#	stdout or to a file.

	###--- Public Methods ---###

	def __init__ (self):
		# Purpose: constructor
		# Returns: nothing
		# Assumes: nothing
		# Effects: initializes 'self', including starting the
		#	timer (for time-stamping)
		# Throws: nothing

		# ordered list of tuples, each being: (time-stamp, message)
		self.messages = []

		# float; number of seconds since the epoch.  This is the
		# starting time for the timer.
		self.initialTime = None

		# log an initial message to get the timer started
		self.stamp ('Profiler initialized')

		return

	def stamp (self,
		message		# string; the message to be time-stamped and
				# ...logged
		):
		# Purpose: time-stamp the given 'message' and add it to the
		#	list of them in 'self'
		# Returns: nothing
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing
		# Notes: starts the internal timer, if it is not already
		#	started

		# if the timer has not been started, then start it by noting
		# the current time as the initial time

		if self.initialTime == None:
			self.initialTime = time.time()
			myTime = 0.0
		else:
			# otherwise, just get the time elapsed so far

			myTime = time.time() - self.initialTime

		self.messages.append ( (myTime, message) )
		return

	def getEntries (self):
		# Purpose: get a list of tuples, one per message logged so
		#	far.
		# Returns: list of tuples, each being:
		#	(float time-stamp, string message)
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing

		return copy.deepcopy (self.messages)

	def write (self,
		filename = None		# optional string; filename to which
					# ...to write the table of profiling
					# ...information
		):
		# Purpose: write a table of profiling information, either to
		#	a file (if 'filename' is specified) or to stdout (if
		#	'filename' is None)
		# Returns: nothing
		# Assumes: we can write to 'filename' if it is specified
		# Effects: writes to stdout or to the file system
		# Throws: IOError if we cannot write to the given 'filename'

		# get a file pointer (fp) to which to write

		if filename == None:
			fp = sys.stdout
		else:
			fp = open (filename, 'w')

		# float; total time measured by the Profiler.  This is the
		# time of the last logged message.
		self.totalTime = 0.0
		if self.messages:
			self.totalTime = self.messages[-1][0]

		# float; time of the last entry we printed
		lastTime = 0.0

		# output the table heading
		self.__printHeading(fp)

		for (timestamp, message) in self.messages:

			# float; time elapsed between the last message printed
			# and the current message
			diffTime = timestamp - lastTime

			# output one or more lines for the current 'message'
			self.__printEntry(fp, timestamp, diffTime, message)

			# remember this time for handling the next message
			lastTime = timestamp

		# if we wrote to a file rather than stdout, then we need to
		# close the file pointer

		if filename != None:
			fp.close()
		return

	###--- Private Methods ---###

	def __printHeading (self,
		fp			# file pointer to which to write
		):
		# Purpose: write the heading for the output table
		# Returns: nothing
		# Assumes: nothing
		# Effects: writes to 'fp'
		# Throws: nothing
		# Notes: Writes a five column table header followed by a line
		#	with hyphens to separate the header row from the data
		#	rows.  The five columns are for elapsed time from the
		#	beginning up to this entry, elapsed time between this
		#	entry and the previous one, percentage of the total
		#	time represented by this elapsed time, a histogram
		#	which prints one '*' for each 5% in this row, and the
		#	actual message string itself.

		line = '%7s %7s %7s %-20s%s\n' % (
			'total',
			'elapsed',
			'percent',
			'histogram',
			'message'
			)
		fp.write(line)
		fp.write(regsub.gsub ('[^ \n]', '-', line))
		return

	def __printEntry (self,
		fp,		# file pointer to which to write
		timestamp,	# float; timestamp for this entry
		diffTime,	# float; amount of time elapsed between this
				# ...entry and the previous one
		message		# string; the message logged
		):
		# Purpose: writes one or more lines needed to represent this
		#	logged entry in the output.
		# Returns: nothing
		# Assumes: nothing
		# Effects: writes to 'fp'
		# Throws: nothing
		# Notes: If 'message' is longer than can be represented in the
		#	rightmost column, then we split it into multiple lines
		#	with those after the first line printed with the first
		#	four columns blank.

		# determine what percentage of the total time was spent for
		# this entry

		if self.totalTime > 0.0:
			percent = diffTime * 100.0 / self.totalTime
		else:
			percent = 0.0

		# integer number of '*' needed, with one per 5% in this entry
		starCount = int(percent / 5.0)

		# actual string of '*' to be printed as the histogram field
		stars = '*' * starCount

		# list of strings; break up long messages as needed into
		# multiple lines
		messageLines = self.__splitMessage(message)

		# print the first line of the message with all five columns
		# filled in

		line = '%7.2f %7.2f %7.2f %-20s%s\n' % (
			timestamp,
			diffTime,
			percent,
			stars,
			messageLines[0]
			)
		fp.write(line)

		# print any following lines with the first four columns blank
		# and only the last filled in

		for messagePart in messageLines[1:]:
			fp.write ('%s%s\n' % (' ' * 44, messagePart) )
		return

	def __splitMessage (self,
		message			# string; the message to be split up
		):
		# Purpose: if 'message' is too long to print in the right
		#	column of one line, then break it up into multiple
		#	strings
		# Returns: list of strings, each of which will fit in the
		#	right-hand column
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing

		# list of strings, each of which is part of 'message'
		list = []

		# integer length of 'message'
		msgLen = len(message)

		# integer; where to start copying the next substring of
		# 'message'
		start = 0

		# integer; how many character will fit in the right column
		chunkSize = 34

		# continue as long as there's more of the 'message' string to
		# be handled
		while start < msgLen:

			# integer; point at which to stop copying this chunk
			# from 'message'
			end = start + chunkSize

			# extract this substring from 'message' and add it to
			# the list

			list.append (message[start:end])
			start = end

		# if the list was empty, then just add an empty string to it

		if len(list) == 0:
			list.append ('')
		return list

###--- Main Program ---###

# If run as a program (instead of being imported as a module), this code will
# run to produce some sample output.

if __name__ == '__main__':

	profiler = Profiler()
	
	list = []
	dict = {}
	maxCount = 10000

	i = 0
	while i < maxCount:
		list.append (i)
		i = i + 1
	profiler.stamp ('created list with %d entries' % maxCount)

	i = 0
	while i < maxCount:
		dict[i] = 1
		i = i + 1
	profiler.stamp ('created dict with %d entries' % maxCount)

	i = 0
	while i < maxCount:
		if i in list:
			pass
		i = i + 1
	profiler.stamp ('tested for %d items in list' % maxCount)

	i = 0
	while i < maxCount:
		if dict.has_key(i):
			pass
		i = i + 1
	profiler.stamp ('tested for %d items in dict' % maxCount)

	profiler.write()
