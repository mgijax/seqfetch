#!/usr/local/bin/python

# Name: MouseMineClient.py
# Purpose: This file provides a python class for interacting with a MouseMine instance
#    to retrieve sequence for canonical and strain-specific genes.  Adding flank to the
#    sequences is not yet supported, as it will require a MouseMine change.
# Notes: This file can be executed as a script for testing purposes; specify IDs as
#    command-line parameters.

import sys
import types
import urllib

class MouseMineClient:
    def __init__ (self, mouseMineUrl):
        # create the client and initialize it with the base URL for MouseMine
        # (eg- 'http://bhmgimm-dev:8080/' or 'http://www.mousemine.org/')

        self.mouseMineUrl = mouseMineUrl
        if not self.mouseMineUrl.endswith('/'):
            self.mouseMineUrl = self.mouseMineUrl + '/'
        return
    
    def getFasta (self, genes):
        # Purpose: retrieve the FASTA for the specified 'genes'.
        # Returns: (string of FASTA sequences, string of failed IDs)
        # Notes: 'genes' can be either a string with a single ID or can be a list of IDs.
        #    Each ID can be either a strain gene ID or can be a canonical gene ID.  If it is
        #    a canonical gene ID, then we convert it to the corresponding C57BL/6J strain
        #    gene ID and request that sequence.
        
        if type(genes) == types.StringType:
            genes = [ genes ]
            
        mapping = self._getCanonicalGeneMapping(genes)
        request = self._buildRequest(genes, mapping)
        return self._processResults(self._convertIDs(self._getResponse(request), mapping), genes)
        
    def _getCanonicalGeneMapping(self, genes):
        # Purpose: get a dictionary that maps from any canonical gene IDs in 'genes' to their
        #    C57BL/6J strain gene equivalents
        
        map = {}
        for gene in genes:
            if gene.startswith('MGI:'):
                parts = gene.split(':')
                map[gene] = '%s_C57BL6J_%s' % (parts[0], parts[1])
        return map

    def _buildRequest(self, genes, mapping):
        # Purpose: Build the XML query string for MouseMine given the list of genes.
        
        request = [
            '<query name="" model="genomic" view="Gene.id">',
            '<constraint path="Gene.primaryIdentifier" op="ONE OF">',
            ]
        
        # add one value for each requested gene, mapping any that have an alternate ID
        for gene in genes:
            if gene in mapping:
                request.append('<value>%s</value>' % mapping[gene])
            else:
                request.append('<value>%s</value>' % gene)

        request.append('</constraint>')
        request.append('</query>')

        return '\n'.join(request)
    
    def _getResponse(self, request):
        # Purpose: Send the given XML 'request' to MouseMine and return its response.

        view = 'Gene.primaryIdentifier'
        url = self.mouseMineUrl + 'mousemine/service/query/results/fasta'
        args = { 'query' : request, 'view' : view }

        fd = urllib.urlopen(url, urllib.urlencode(args))
        lines = fd.readlines()
        fd.close()

        return ''.join(lines)
    
    def _convertIDs(self, response, mapping):
        # Purpose: Go through the 'response' from MouseMine and convert any IDs that we mapped
        #    back to their original versions.
        
        for (canonical, strainSpecific) in mapping.items():
            response = response.replace(strainSpecific, canonical)
        
        return response
    
    def _processResults(self, response, genes):
        # Purpose: Sort the results in 'response' so that they are in the same order as the IDs
        #   used in the request (in 'genes').  If any gene IDs are not in the result set, add
        #   a message about them.
        # Returns: (string of FASTA sequences found and ordered by IDs in 'genes',
        #   string with any IDs from 'genes' for which we found no sequence)
        
        responseMap = {}                # maps ID to FASTA data
        sequences = response.split('>')
        for seq in sequences[1:]:
            accID = seq[:seq.find(' ')]
            responseMap[accID] = '>' + seq
            
        out = []
        failed = []

        for geneID in genes:
            if geneID in responseMap:
                out.append(responseMap[geneID])
            else:
                failed.append(geneID)

        if failed:
            return ''.join(out), 'MouseMine could not find: %s' % ', '.join(failed)
        return ''.join(out), ''
    
# execute this file as a script to allow interactive testing
if __name__ == '__main__':
    if len(sys.argv) > 1:
        mmc = MouseMineClient('http://bhmgimm-dev:8080/')
        print mmc.getFasta(sys.argv[1:])