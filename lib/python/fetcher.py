# Name: fetcher.py
# Purpose: Retrieves various types of sequences from online resources, allowing us to stop
#    maintaining those data sets locally.

import sys
from urllib.request import urlopen
from urllib.parse import urlencode
import json

# default values for build and strain (can override)
genomeBuild = 'GRCm38.p6'
mouseStrain = 'C57BL/6J'

###--- functions ---###

def setGenomeBuild(build):
    # set the module's genome build
    global genomeBuild
    genomeBuild = build
    return

def setMouseStrain(strain):
    # set the module's mouse strain
    global mouseStrain
    mouseStrain = strain
    return

###--- classes ---###

# Base class for fetching sequences, not to be instantiated directly.
class SequenceFetcher:
    # maps a base to its complement
    BASE_COMPLEMENT = { 
      'a' : 't',
      't' : 'a',
      'c' : 'g',
      'g' : 'c',
      'n' : 'n',
      
      'A' : 'T',
      'T' : 'A',
      'C' : 'G',
      'G' : 'C',
      'N' : 'N' 
    }

    def complement (self, dna) :
        # Returns the complement of sequence 'dna'.
        # Assumes all characters in 'dna' appear as keys in BASE_COMPLEMENT.
        return ''.join([self.BASE_COMPLEMENT.get(b, b) for b in dna])

    def reverseComplement (self, dna) :
        # Returns the reverse complement of sequence 'dna'.  That is, it complements the 'dna'
        # string and then reverses it for the minus (-) strand.
        # Assumes all characters in 'dna' appear as keys in BASE_COMPLEMENT.
        return self.complement(dna)[::-1]

    def chunkString (self, s, n) :
        # Breaks string 's' up into lines of up to 'n' characters each.
        # Returns the resulting list of strings.
        return [s[i:i+n] for i in range(0, len(s), n)]

    def fetchById(self, id):
        # Returns the sequence corresponding to the given seq 'id'.
        return self._fetch(self.BASEURL % id)

    def _fetch (self, url, args = None) :
        # Read from the given 'url' (and passing along any extra 'args').
        # Returns the string that is read.
        if args:
            fd = urlopen(url, urlencode(args).encode('ascii'))
        else:
            fd = urlopen(url)
        s = fd.read()
        s = s.decode('utf-8')
        fd.close()
        return s

# Is a SequenceFetcher for reading from the UniProt resource.
class UniprotFetcher (SequenceFetcher) :
    BASEURL="https://www.uniprot.org/uniprot/%s.fasta"

# Is a SequenceFetcher for reading from the Entrez resource at NCBI.
class EntrezFetcher (SequenceFetcher) :
    BASEURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=fasta&retmode=text"

# Is a SequenceFetcher for reading from the Ensembl resource.
class EnsemblFetcher (SequenceFetcher) :
    BASEURL = "http://rest.ensembl.org/sequence/id/%s?content-type=text/x-fasta"

# Is a SequenceFetcher for reading from the MouseMine resource at MGI.
class MouseMineFetcher (SequenceFetcher) :
    BASEURL = "http://www.mousemine.org/mousemine/service/"

    def fetchById(self, id):
        # Returns the sequence corresponding to the given seq 'id'.
        url = self.BASEURL + "query/results/fasta"
        args = {
            'query' : '''
                <query model="genomic" view="SequenceFeature.primaryIdentifier" >
                <constraint path="SequenceFeature.primaryIdentifier" op="=" value="%s"/>
                </query>
                ''' % id,
            'view' : 'SequenceFeature.primaryIdentifier'
        }
        return self._fetch(url, args)

    def fetchByCoordinates (self, build, strain, chrom, start, end, strand, flank = 0) :
        # Returns a slice of the genomic sequence corresponding to the given input parameters.
        # (genome build, strain, chromosome, start coordinate, end coordinate, strand, and the amount of flank to include)
        url = self.BASEURL + "sequence"
        args = {
            'start' : max(0, start - flank - 1),
            'end' : end + flank,
            'query' : '''
                <query model="genomic" view="Chromosome.sequence.residues">
                  <constraint path="Chromosome.strain.name" op="=" value="%s" />
                  <constraint path="Chromosome.primaryIdentifier" op="=" value="%s" />
                  </query>
                ''' % (strain, chrom)
        }

        # Read from the given 'url' and convert the resulting JSON string into a Python dictionary (and associated structures).
        result = json.loads(self._fetch(url, args))

        # Extract the sequence, reverse complement it if the minus strand was requested, and break it into 60-character chunks.
        seq = result['features'][0]['seq']
        if strand == "-":
            seq = self.reverseComplement(seq)
        seq = self.chunkString(seq, 60)

        # Add a header line to make complete the FASTA format, and return the result.
        hdr = ">dna/%s/chr%s:%d..%d(%s)" %(build, chrom, start, end, strand)
        return hdr + '\n' + '\n'.join(seq) + '\n'


# Maps from a sequence database type to the class that should be used to fetch its sequences.
type2class = {
    "swissprot" : UniprotFetcher,
    "trembl"    : UniprotFetcher,
    "sptrembl"  : UniprotFetcher,
    "genbank"   : EntrezFetcher,
    "refseq"    : EntrezFetcher,
    "ensembl_mus_cdna" : EnsemblFetcher,
    "ensembl_mus_prot" : EnsemblFetcher,
    "straingene" : MouseMineFetcher,
    "mousegenome" : MouseMineFetcher,
}

def fetch (arg) :
    # Fetch the sequence corresponding to the given sequence identification string (in 'arg').
    # Throws an Exception in case of any failure
    
    db,id,chr,start,end,strand,flank = arg.split("!")
    cls = type2class[db]
    fetcher = cls()
    
    # If no start coordinate, assume we will fetch by ID.
    if start == '':
        return fetcher.fetchById(id)

    else:
        return fetcher.fetchByCoordinates(
            genomeBuild, mouseStrain, chr, int(start), int(end), strand, (int(flank) if flank else 0))

def _test_ () :
    # Run automated tests using a set of pre-defined sequence identification strings, writing
    # to stdout.
    
    testcases = [
        "mousegenome!test!10!100015630!100015640!+!5",
        "mousegenome!test!10!100015630!100015640!-!5",
        "mousegenome!ENSMUSG00000015452!17!34597862!34600936!+!",
        "mousegenome!ENSMUSG00000015452!17!34597862!34600936!+!1000",
        "ensembl_mus_cdna!ENSMUST00000105283!!!!!",
        "ensembl_mus_prot!ENSMUSP00000100920!!!!!",
        "refseq!XM_006513314!!!!!",
        "genbank!AK134301!!!!!",
        "swissprot!P20826!!!!!",
        "trembl!A0A1W2P773!!!!!",
        "straingene!MGP_LPJ_G0017604!!!!!",
        "straingene!MGP_C3HHeJ_G0017455!!!!!",
    ]
    print("===============Automated Tests==============")
    for arg in testcases:
        print("\n============================================")
        print(arg)
        print(fetch(arg))

# Invoke as a script with no arguments for automated tests or with a series of sequence
# identification strings to retrieve them and report to stdout.
if __name__ == "__main__":
     if len(sys.argv) > 1:
         for arg in sys.argv[1:]:
             print(fetch(arg))
     else:
         _test_()

