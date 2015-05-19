#!/usr/bin/python
# 
# 
# 
# 
# Kim Brugger (30 Apr 2015), contact: kim@brugger.dk

import sys
import pprint
pp = pprint.PrettyPrinter(indent=4)
import re
import scipy
from scipy import stats
import getopt

sys.path.append("/software/lib/python2.7/site-packages/pysam-0.7.5-py2.7-linux-x86_64.egg")
import pysam


#
# reverse a DNA string, even with IUPAC codes
#
# Kim Brugger (30 Apr 2015), contact: kbr@brugger.dk
def revDNA( string ):

    string = string.upper()
    
    rev_bases = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  '-': '-', 'N': 'N',
                  'W': 'S', 'S': 'W', 'M': 'K', 'K': 'M',
                  'R': 'Y', 'Y': 'R', 'B': 'B', 'D': 'D',
                  'H': 'H', 'V': 'V'}

    rev_str = len(string)*[None]
    for i in range(0, len(string)):
        rev_str[ len(string) - i - 1] =  rev_bases[ string[ i ]]

    return "".join( rev_str )


#
# creates an IUPAC based on the input of 2-4 seqs input
def IUPAC_consensus( bases ):

    # This bit is shit, need to fix it later...
    fb = ""
    sb = "" 
    tb = "" 
    Fb = ""

    if ( len(bases) >= 1):
        fb = bases[0]
    if ( len(bases) >= 2):
        sb = bases[1] 
    if ( len(bases) >= 3):
        tb = bases[2]
    if ( len(bases) >= 4):
        Fb = bases[3]

    max_length =  max(len(fb), len(sb), len(tb), len(Fb))

    IUPAC_Codes = {'A'   : 'A',
                   'C'   : 'C',
                   'G'   : 'G',
                   'T'   : 'T',
                   'AT'  : 'W',
                   'CG'  : 'S',
                   'AC'  : 'M',
                   'GT'  : 'K',
                   'AG'  : 'R',
                   'CT'  : 'Y',
                   'CGT' : 'B',
                   'AGT' : 'D', 
                   'ACT' : 'H',
                   'ACG' : 'V',
                   'ACGT': 'N',

                   '-'    : '-',
                   'A-'   : 'A',
                   'C-'   : 'C',
                   'G-'   : 'G',
                   'T-'   : 'T',
                   'AT-'  : 'W',
                   'CG-'  : 'S',
                   'AC-'  : 'M',
                   'GT-'  : 'K',
                   'AG-'  : 'R',
                   'CT-'  : 'Y',
                   'CGT-' : 'B',
                   'AGT-' : 'D', 
                   'ACT-' : 'H',
                   'ACG-' : 'V',
                   'ACGT-': 'N'

                   }

    IUPAC_consensus = ""

    for i in range(0, max_length):
        bases = [];

        if (len(fb) > i):
            bases.append( fb[i] )

        if (len(sb) > i):
            bases.append( sb[i] )

        if (len(tb) > i):
            bases.append( tb[i] )

        if (len(Fb) > i):
            bases.append( Fb[i] )

        fst = ''.join( sorted(set(bases)))

        if ( fst != ""):
            IUPAC_consensus += IUPAC_Codes[ fst ]

    return(IUPAC_consensus)

#
# Goes from a 3bp seq to a codon. Surprise!!!!
# 
# Kim Brugger (30 Apr 2015), contact: kbr@brugger.dk
def codon2AA(codon):

    codon = codon.upper()

    codon2AA = { 'TTT': 'F', 'TTC': 'F',
                 'TTA': 'L', 'TTG': 'L',
                 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                 'TAT': 'Y', 'TAC': 'Y', 
                 'TAA': '*', 'TAG': '*', 'TGA': '*',
                 'TGT': 'C', 'TGC': 'C',
                 'TGG': 'W',
                 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                 'CAT': 'H', 'CAC': 'H',
                 'CAA': 'Q', 'CAG': 'Q',
                 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
                 'ATG': 'M',
                 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                 'AAT': 'N', 'AAC': 'N',
                 'AAA': 'K', 'AAG': 'K',
                 'AGT': 'S', 'AGC': 'S',
                 'AGA': 'R', 'AGG': 'R',
                 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                 'GAT': 'D', 'GAC': 'D',
                 'GAA': 'E', 'GAG': 'E',
                 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G' }
                 
    if (len(codon) < 3 or re.search('-', codon)):
        del_len = codon.count('-');
        return "d %d bp" % del_len;

    if (len(codon) > 3):
        return "i %d bp" %( len(codon) - 3)

    if ( codon not in codon2AA ):
#        print "Cannot tranlate " + codon + " codon to an AA";
        return "X"

    return codon2AA[ codon ] 


#
# dict for storing the stats regarding the bases we are currently looking at
#
# Kim Brugger (30 Apr 2015), contact: kbr@brugger.dk
def update_base_counts(base_counts, base, qual, strand, pos, read_nr):

    if ( base not in base_counts):
        base_counts[ base ] = dict()
        base_counts[ base ][ 'count' ] = 0
        base_counts[ base ][ 'quals' ] = 0
        base_counts[ base ][ 'read1' ] = 0
        base_counts[ base ][ 'read2' ] = 0
        base_counts[ base ][ 'starts' ] = []
        # plus and minus strand
        base_counts[ base ][ 0 ] = dict()
        base_counts[ base ][ 1 ] = dict()

    base_counts[ base ][ 'count' ] += 1
    base_counts[ base ][ 'quals' ] += qual
    base_counts[ base ][ 'starts' ].append( pos )

    if ( read_nr == 1 ):
        base_counts[ base ][ 'read1' ] += 1
    else:
        base_counts[ base ][ 'read2' ] += 1

    if (pos not in base_counts[ base ][ strand ]):
        base_counts[ base ][ strand ][ pos ] = 1
    else:
        base_counts[ base ][ strand ][ pos ] += 1

    return base_counts;


#
# dict for storing the stats regarding the codons we are currently looking at
#
# Kim Brugger (30 Apr 2015), contact: kbr@brugger.dk
def update_AA_counts(AA_counts, codon, qual, strand, pos, read_nr):

    AA = codon2AA( codon )


    if (AA not in AA_counts):
        AA_counts[ AA ] = dict()
        AA_counts[ AA ][ 'count' ] = 0
        AA_counts[ AA ][ 'quals' ] = 0
        AA_counts[ AA ][ 'read1' ] = 0
        AA_counts[ AA ][ 'read2' ] = 0
        # plus and minus strand
        AA_counts[ AA ][ 0 ] = dict()
        AA_counts[ AA ][ 1 ] = dict()


    AA_counts[ AA ][ 'count' ] += 1
    AA_counts[ AA ][ 'quals' ] += qual

    if ( read_nr == 1 ):
        AA_counts[ AA ][ 'read1' ] += 1
    else:
        AA_counts[ AA ][ 'read2' ] += 1

    if (pos not in AA_counts[ AA ][ strand ]):
        AA_counts[ AA ][ strand ][ pos ]  = 1
    else:
        AA_counts[ AA ][ strand ][ pos ] += 1

    return AA_counts;




#
# Sets starts and stop positions if provided by the user, otherwise do the whole sequence
#
# Kim Brugger (30 Apr 2015), contact: kbr@brugger.dk
def set_start_n_stop():
    global start, enda

    if ( len(sys.argv) >= 4 ) :
        start = int(sys.argv[3])


    if ( len(sys.argv) >= 5 ) :
        end = int(sys.argv[4])


    if ( start > -1 and end > -1 and (end < start)):
        minus_strand = 1;
        (start, end) = (end, start)



def set_gene_list():

    global fasta_ref, gene_list

    if ( re.search('C06', bamfile)):
        gene_list = [['Pol', 2253, 2546, 5],
                     ['RT', 2550, 4227, 0, 320],
                     ['Int',4233, 5100, 1]]

        fasta_ref = pysam.Fastafile("/refs/HIV/K03455.fasta");
        reference = "K03455";

    elif ( re.search(r'C11', bamfile)):
    # The gene range needs to include the stop codon on the minus strand.
        gene_list = [['TK/UL23', 48000, 46873],
                     ['Pol/UL30', 63265, 66987]]

        fasta_ref = pysam.Fastafile("/refs/HSV/Z86099.fasta");
        reference = "Z86099";
        
    elif ( re.search(r'C19', bamfile)):
        gene_list = [['Kinase', 140484, 142405], 
                     ['Pol', 80631,76906]]
        fasta_ref = pysam.Fastafile('/refs/CMV/CMV_AD169.fasta')
        reference = "CMV_AD169"
    else:
        gene_list = [['Pol', 2253, 2546, 5],
                     ['RT', 2550, 4227, 0, 320],
                     ['Int',4233, 5100, 1]]

        fasta_ref = pysam.Fastafile("/refs/HIV/K03455.fasta");
        reference = "K03455";


def average_base_qual(base_quals):

    summed_quals = 0

    for base_qual in base_quals:
        summed_quals += ord( base_qual ) - 33

    base_qual = summed_quals/len( base_quals )

    return base_qual


# ---------------- MAIN LOOP ------------------------
bamfile   = sys.argv[1]
fasta_ref = sys.argv[2]

bam       = pysam.Samfile( bamfile, "rb" )


start        = -1
end          = -1
minus_strand = 0
set_start_n_stop()
#print str(start) + " " + str(end)

gene_list = []
set_gene_list()
#pp.pprint( gene_list )

MIN_MAPQ          =   20
MIN_BASEQ         =   30
MIN_COVERAGE      = 1000
MIN_MAP_LEN       =   80
MIN_BP_FROM_END   =    5
MIN_ALLELE_PERC   =   20
AA_MIN_PERC       =    1



Stanford_format = dict()
fasta_consensus = []
AA_changes      = []

ORF_start = 2550
ORF_end   = 4227
gene_name = 'Kinase'
Stanford_format[gene_name] = [];


ORF_start = 140484
ORF_end   = 142405

def find_consensus_base( base_counts ):

    base_freqs = dict()

    total_bases = 0
    for base in base_counts:
        total_bases += base_counts[ base ]['count']

    for base in base_counts:
        base_freq = base_counts[ base ]['count']*100.00/total_bases;
        if ( base_freq < MIN_ALLELE_PERC ):
            continue

        base_freqs[ base ] = base_freq

#    pp.pprint( base_freqs )

    if ( total_bases < MIN_COVERAGE ):
        return 'N'


    if (len( base_freqs.keys()) >= 1):

        bases_by_freq = sorted(base_freqs, key=base_freqs.get, reverse=True)

#        pp.pprint( base_counts )

        for i in range(1, len( bases_by_freq)):

            major_base = bases_by_freq[ 0 ]
            minor_base = bases_by_freq[ i ]

            major_freq = base_freqs[ major_base ]
            minor_freq = base_freqs[ minor_base ]

            major_starts_fwd = len(base_counts[ major_base ][ 0 ].keys())
            minor_starts_fwd = len(base_counts[ minor_base ][ 0 ].keys())

            major_starts_rev = len(base_counts[ major_base ][ 1 ].keys())
            minor_starts_rev = len(base_counts[ minor_base ][ 1 ].keys())

            major_starts =  major_starts_fwd + major_starts_rev
            minor_starts =  minor_starts_fwd + minor_starts_rev

            major_read1 = base_counts[ major_base ]['read1'] 
            major_read2 = base_counts[ major_base ]['read2']

            minor_read1 = base_counts[ minor_base ]['read1'] 
            minor_read2 = base_counts[ minor_base ]['read2']


            print "\t".join([str(pile.pos),
                             major_base,
                             minor_base,
                             str(base_counts[ major_base ][ 'count' ]),
                             str(base_counts[ minor_base ][ 'count' ]),
                             "%.2f" % major_freq,
                             "%.2f" % minor_freq,
                             str(major_starts),
                             str(major_starts_fwd),
                             str(major_starts_rev),
                             str(minor_starts),
                             str(minor_starts_fwd),
                             str(minor_starts_rev),
                             
                             str(major_read1),
                             str(major_read2),
                             
                             str(minor_read1),
                             str(minor_read2)])
            continue
                             

            if ( minor_freq < 3):
                print "%d -> %d [%s %s] [%f %f]" % (0, i +1 , major_base, minor_base, major_freq, minor_freq)
 #               print "%d %d <--> %d" % (pile.pos, base_counts[ major_base ]['read1'], base_counts[ major_base ]['read2'])
                print "STRANDS [%d] %d  %d" % (pile.pos, base_counts[ minor_base ]['read1'], base_counts[ minor_base ]['read2'])

                print "\t".join([str(len(base_counts[ major_base ][ 0 ].keys())),
                                 str(len(base_counts[ major_base ][ 1 ].keys())),
                                 str(len(base_counts[ minor_base ][ 0 ].keys())),
                                 str(len(base_counts[ minor_base ][ 1 ].keys()))])

            major_starts =  (base_counts[ major_base ][ 0 ].keys() + base_counts[ major_base ][ 1 ].keys())
            minor_starts =  (base_counts[ minor_base ][ 0 ].keys() + base_counts[ minor_base ][ 0 ].keys())


            start_bias_zvalue, start_bias_pvalue_fwd = stats.ranksums(major_starts, minor_starts)

#            print abs(start_bias_zvalue)
#            sleep( 1 )


            start_bias_zvalue_fwd, start_bias_pvalue_fwd = stats.ranksums(major_starts_fwd, minor_starts_fwd)


            major_starts_rev =  (base_counts[ major_base ][ 1 ].keys())
            minor_starts_rev =  (base_counts[ minor_base ][ 1 ].keys())
            start_bias_zvalue_rev, start_bias_pvalue_rev = stats.ranksums(major_starts_rev, minor_starts_rev)


            if ( start_bias_pvalue_fwd < 0.05 or start_bias_pvalue_fwd < 0.05):
                print str(start_bias_pvalue_fwd) + " " + str(start_bias_pvalue_rev)
                print "%s (%f) vs %s (%f) fails, should it drop %s?" % ( major_base, major_freq, minor_base, minor_freq, minor_base)                
#                del base_freqs[ minor_base]
#                exit()

            # Use fishers test to see if there is a strand bias for the SNP
            oddsratio, strand_bias_pvalue = stats.fisher_exact([[len(base_counts[ major_base ][ 0 ].keys()), 
                                                                 len(base_counts[ major_base ][ 1 ].keys())], 
                                                                [len(base_counts[ minor_base ][ 0 ].keys()), 
                                                                 len(base_counts[ minor_base ][ 1 ].keys())]])
            if ( strand_bias_pvalue < 0.05 ):
#                pp.pprint( base_counts)
                print "\t".join([str(len(base_counts[ major_base ][ 0 ].keys())),
                                 str(len(base_counts[ major_base ][ 1 ].keys())),
                                 str(len(base_counts[ minor_base ][ 0 ].keys())),
                                 str(len(base_counts[ minor_base ][ 1 ].keys()))])

#                print "%d -> %d [%s %s] [%f %f]" % (0, i +1 , major_base, minor_base, major_freq, minor_freq)
#                print "strand bias (pvalue): %.2f (> 0.05 no bias)" % strand_bias_pvalue
#                exit()


    

#    print pile.pos
    consensus_base = IUPAC_consensus( (base_freqs.keys()))

    return consensus_base



def find_significant_AAs( AA_counts, genome_pos, ref_AA, AA_number ):

    AA_freqs = dict()

    total_AAs = 0
    for AA in AA_counts:
        total_AAs += AA_counts[ AA ]['count']

    if ( total_AAs < MIN_COVERAGE ):
        return 'x'
    
#    pp.pprint( AA_counts )    

    alternative_AAs = 0

    for AA in AA_counts:
        if ( ref_AA == AA ):
            continue

        AA_freq = AA_counts[ AA ]['count']*100.00/total_AAs;
        if ( AA_freq < AA_MIN_PERC ):
            continue

        AA_freqs[ AA ] = AA_freq
        alternative_AAs += 1

    if ( alternative_AAs == 0 ):
        return

    AAs_by_freq = sorted(AA_freqs, key=AA_freqs.get, reverse=True)
    AA_line = []
    AA_line.append(str(genome_pos));
    AA_line.append(ref_AA+str(AA_number))
    for AA in AAs_by_freq:
        AA_line.append( "%s:%.2f%%" % (AA, AA_freqs[ AA ]))
        Stanford_format[gene_name].append("%s%d%s" % (ref_AA, AA_number, AA))
        
    AA_changes.append("\t".join( AA_line ))
#    pp.pprint( AA_counts )
#    exit()
    

if ( 0 ) :
    print "\t".join(["pile.pos",
                 "major_base",
                 "minor_base",
                 "major count",
                 "minor count",

                 "major_freq",
                 "minor_freq",
                 "major_starts",
                 "major_starts_fwd",
                 "major_starts_rev",
                 "minor_starts",
                 "minor_starts_fwd",
                 "minor_starts_rev",
                         
                 "major_read1",
                 "major_read2",
                 
                 "minor_read1",
                 "minor_read2"])


deletion_skipping = 0

#for pile in bam.pileup():
#for pile in bam.pileup('K03455', 2670, 4227,max_depth=1000):
#for pile in bam.pileup('K03455', 2263, 2265,max_depth=10000000):
max_depth = 0
for pile in bam.pileup('CMV_AD169', 142272, 142291,max_depth=10000):



    if ( start > -1 and end > -1) :
        if (pile.pos < start or pile.pos > end):
            continue

    passed_bases = 0

    base_counts  = dict()
    AA_counts = dict()
    insertion_count = 0
    deletion_count  = 0

    # We are in a coding frame...
    if ( (ORF_start - 1  - pile.pos )%3 == 0):
            
        ref_codon = fasta_ref.fetch(str(bam.getrname(pile.tid)), pile.pos, pile.pos+3 )
        ref_AA    = codon2AA( ref_codon )


    if (max_depth < pile.n ):
        max_depth = pile.n


    for read in pile.pileups:

        if (read.alignment.is_unmapped or read.alignment.is_duplicate or read.is_del):
                continue

        if (read.alignment.mapq <= MIN_MAPQ ):
            continue

        if (read.qpos < MIN_BP_FROM_END or read.qpos + MIN_BP_FROM_END > read.alignment.qend):
            continue

        if ( read.alignment.alen < MIN_MAP_LEN ):
            continue

        read_nr = 1
        if ( read.alignment.is_read2):
            read_nr = 2

        avg_base_qual = 0
        alt           = "N"

        if ( read.indel > 0):
            alt = read.alignment.seq[ read.qpos:read.qpos+ read.indel + 1]
            avg_base_qual = average_base_qual(read.alignment.qual[ read.qpos:read.qpos+ read.indel + 1]);

        elif ( read.indel < 0):
            alt = fasta_ref.fetch(str(bam.getrname(pile.tid)), pile.pos, pile.pos+abs(read.indel)+1 )

            alt = read.alignment.seq[ read.qpos ] + "-" * abs(read.indel)

            avg_base_qual = ord( read.alignment.qual[ read.qpos] ) - 33

        else:
            avg_base_qual = ord( read.alignment.qual[ read.qpos] ) - 33
            alt = read.alignment.seq[ read.qpos];

        if ( avg_base_qual > MIN_BASEQ):
            base_counts = update_base_counts(base_counts, alt, 
                                             avg_base_qual, read.alignment.is_reverse, read.alignment.aend, read_nr)

        # We are in a coding frame...
        if ( (ORF_start - 1  - pile.pos )%3 == 0):
            
            if ( read.qpos + 3 > read.alignment.qend):
                print "going across the end ... %d %d" % ( read.qpos, read.alignment.qend )
                continue

            codon          = "NNN"
            avg_codon_qual = 0
 
            if ( read.indel > 0):
                insertion_count += 1
                avg_codon_qual = average_base_qual(read.alignment.qual[ read.qpos:read.qpos +  3]);
                codon = read.alignment.seq[ read.qpos:read.qpos + 3 ] + '+'

            elif ( read.indel < 0):
                avg_codon_qual = average_base_qual(read.alignment.qual[ read.qpos:read.qpos +  3]);
                deletion_count += 1
                codon = read.alignment.seq[ read.qpos] + "-" * abs(read.indel) + read.alignment.seq[ read.qpos + 1: read.qpos + 2]
#                if (len(codon)>3):
#                    codon = codon[0:3]

            else:
                codon = read.alignment.seq[ read.qpos:read.qpos + 3];
                avg_codon_qual = average_base_qual(read.alignment.qual[ read.qpos:read.qpos +  3]);

            
            if ( avg_codon_qual > MIN_BASEQ and pile.n > MIN_COVERAGE):
                AA_counts = update_AA_counts(AA_counts, codon, 
                                             avg_codon_qual, read.alignment.is_reverse, read.alignment.aend, read_nr)

        

    if ( (ORF_start - 1  - pile.pos )%3 == 0):

        AA_number    = (1+ (pile.pos + 1 - ORF_start)/3)
        genome_pos      = pile.pos+1;

        find_significant_AAs( AA_counts, genome_pos, ref_AA, AA_number )
#        print ref_AA
#        pp.pprint( codon_counts )
#        exit()


    if ( deletion_skipping > 0 ):
        deletion_skipping -= 1
    else:
        consensus_base = find_consensus_base( base_counts )
        # If the consensus base contains a deletion, set the level of skipping here.
        deletion_skipping =  consensus_base.count('-')

#        print str(pile.pos) + " " + consensus_base+" - (" + str( deletion_skipping ) + ")"

        fasta_consensus.append( consensus_base )




#    if ( pile.pos == 142271):
#        pp.pprint( base_counts )
#        exit()

#    pp.pprint( codon_counts )


print "\n".join( AA_changes )

print " ".join(Stanford_format[ gene_name])


    #continue

#    print consensus_base( base_counts )

    

#if (len(consensus) == 0 or len(consensus) == N_bases):
#    exit()



            
#### get input filename without the path
b = re.search(".*\/?(.*)$", bamfile)
#print b
if ( b.group(1) ):
    bamfilename = b.group(1)
else:
    bamfilename = bamfile



 


#print ">IUPAC_Consensus"


consensus_seq = "".join( fasta_consensus )


if ( minus_strand ):
    consensus_seq = revDNA( consensus_seq )


consensus_seq = re.sub(r'^N+',"", consensus_seq)
consensus_seq = re.sub(r'N+$',"", consensus_seq)

if ( len(consensus_seq) == 0):
    exit()


if ( start > -1 and end > -1):
    print ">%s_Consensus %d-%d" %(bamfilename, start, end)
else:
    print ">%s_Consensus" %(bamfilename)

print consensus_seq


bam.close()

