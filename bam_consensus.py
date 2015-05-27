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


# Various global variables 

FASTA_OUT = 1
CODON_OUT = 2
TAB_OUT   = 4

# These affect the sensitivity and behaviour of  the program
MIN_MAPQ          =    20
MIN_BASEQ         =    30
MIN_COVERAGE      =   400
DOWNSAMPLE_DEPTH  = 10000
MIN_MAP_LEN       =    80
MIN_BP_FROM_END   =     5
MIN_ALLELE_PERC   =    20
AA_MIN_PERC       =     1


#
# reverse a DNA string, even with IUPAC codes and gaps!
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
# creates an IUPAC based on the input of 2-4 seqs input, including gaps
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
                 
    if (re.search('-', codon)):
        del_len = codon.count('-');
        if (del_len > 3):
            del_len = 3

        return "deletion %d bp" % del_len;

    if (len(codon) > 3):
        return "insertion %d bp" %( len(codon) - 3)

    if ( codon not in codon2AA ):
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
# Grow codons base by base to ensure that we capture indels. Damn, this is getting complicated.
#
def grow_codons( codons, read_name, alt, avg_base_qual, strand, pos, read_nr, first_base):

    if ( read_name not in codons ):
        

        #Only if we are on the first base of a codon do we create a new entry
        if ( not first_base ):
            return codons

        codons[ read_name ] = dict()
        codons[ read_name ][ 'codon' ]   = ""
        codons[ read_name ][ 'quals' ]   = 0

        codons[ read_name ][ 'start' ]   = 0
        codons[ read_name ][ 'strand' ]  = 0
        codons[ read_name ][ 'read_nr' ] = 0


    # The two reads overlap at this position so discard thm to remove any bias
    if ( codons[ read_name ][ 'read_nr' ] and codons[ read_name ][ 'read_nr' ] != read_nr ):
        codons.pop( read_name )
        return codons

    codons[ read_name ][ 'codon' ]  +=  alt
    codons[ read_name ][ 'quals' ]  +=  avg_base_qual
    codons[ read_name ][ 'start' ]   = pos
    codons[ read_name ][ 'strand' ]  = strand
    codons[ read_name ][ 'read_nr' ] = read_nr

    return codons


def codons2AA_counts( codons ):

    depth = len(codons.keys())

#    print depth

    AA_counts = dict()

#    pp.pprint( codons )

    for name in codons.keys():

        codon   = codons[ name ][ 'codon' ]
        qual    = codons[ name ][ 'quals' ]
        strand  = codons[ name ][ 'strand' ]
        pos     = codons[ name ][ 'start' ]
        read_nr = codons[ name ][ 'read_nr' ]

        
        # The seq is shorter than the full codon, so we drop it from the table and move on
        if ( len( codons[ name ][ 'codon' ]) < 3):
            codons.pop( name )
        else:

            AA_counts = update_AA_counts(AA_counts, codon, qual, strand, pos, read_nr)


            # Trim away the bases we have used, This is done to ensure that we can 
            # report longer deletions in a nice manner. I am sure there is a massive 
            # computational overhad for this, but it works really nice.
            codons[ name ][ 'codon' ] = codons[ name ][ 'codon' ][3:-1]

            if ( len( codons[ name ][ 'codon' ]) == 0):
                codons.pop(name)

    return (AA_counts, codons)


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
# Finds the consensus base(s)
#
# Kim Brugger (27 May 2015), contact: kbr@brugger.dk
def find_consensus_base( base_counts, pos = 0, ref_base = 'N' ):

    base_freqs = dict()
    all_base_freqs = dict()

    total_bases = 0
    for base in base_counts:
        total_bases += base_counts[ base ]['count']

    for base in base_counts:
        base_freq = base_counts[ base ]['count']*100.00/total_bases;
        all_base_freqs[ base ] = base_freq
        if ( base_freq < MIN_ALLELE_PERC ):
            continue

        base_freqs[ base ] = base_freq

#    pp.pprint( base_freqs )

    if ( total_bases < MIN_COVERAGE ):
        return 'N'

    if (output_format & TAB_OUT and len( all_base_freqs.keys()) >= 1):

        bases_by_freq = sorted(all_base_freqs, key=base_freqs.get, reverse=True)

        if (ref_base == 'N'):
            
            ref_base = bases_by_freq[0]

        ref_read1 = base_counts[ ref_base ]['read1'] or 0
        ref_read2 = base_counts[ ref_base ]['read2']
        ref_starts_rev = len(base_counts[ ref_base ][ 1 ].keys())
        ref_starts_fwd = len(base_counts[ ref_base ][ 0 ].keys())
        ref_freq = all_base_freqs[ ref_base ]
        ref_starts =  ref_starts_fwd + ref_starts_rev

        for i in range(0, len( bases_by_freq)):

            alt_base = bases_by_freq[ i ]

            if ( ref_base == alt_base):
                continue

            alt_freq = all_base_freqs[ alt_base ]


            alt_starts_fwd = len(base_counts[ alt_base ][ 0 ].keys())
            alt_starts_rev = len(base_counts[ alt_base ][ 1 ].keys())
            alt_starts =  alt_starts_fwd + alt_starts_rev
            alt_read1 = base_counts[ alt_base ]['read1'] 
            alt_read2 = base_counts[ alt_base ]['read2']

            global tab_lines
            tab_lines.append( "\t".join([str(pos),
                             ref_base,
                             base,
                             str(base_counts[ ref_base ][ 'count' ]),
                             str(base_counts[ alt_base ][ 'count' ]),
                             "%.2f" % ref_freq,
                             "%.2f" % alt_freq,
                             str(ref_starts),
                             str(ref_starts_fwd),
                             str(ref_starts_rev),
                             str(alt_starts),
                             str(alt_starts_fwd),
                             str(alt_starts_rev),
                             
                             str(ref_read1),
                             str(ref_read2),
                             
                             str(alt_read1),
                             str(alt_read2)]))
            continue

            if ( alt_freq < 3):
                print "%d -> %d [%s %s] [%f %f]" % (0, i +1 , ref_base, alt_base, ref_freq, alt_freq)
 #               print "%d %d <--> %d" % (pile.pos, base_counts[ ref_base ]['read1'], base_counts[ ref_base ]['read2'])
                print "STRANDS [%d] %d  %d" % (pile.pos, base_counts[ alt_base ]['read1'], base_counts[ alt_base ]['read2'])

                print "\t".join([str(len(base_counts[ ref_base ][ 0 ].keys())),
                                 str(len(base_counts[ ref_base ][ 1 ].keys())),
                                 str(len(base_counts[ alt_base ][ 0 ].keys())),
                                 str(len(base_counts[ alt_base ][ 1 ].keys()))])

            ref_starts =  (base_counts[ ref_base ][ 0 ].keys() + base_counts[ ref_base ][ 1 ].keys())
            alt_starts =  (base_counts[ alt_base ][ 0 ].keys() + base_counts[ alt_base ][ 0 ].keys())


            start_bias_zvalue, start_bias_pvalue_fwd = stats.ranksums(ref_starts, alt_starts)
            start_bias_zvalue_fwd, start_bias_pvalue_fwd = stats.ranksums(ref_starts_fwd, alt_starts_fwd)


            ref_starts_rev =  (base_counts[ ref_base ][ 1 ].keys())
            alt_starts_rev =  (base_counts[ alt_base ][ 1 ].keys())
            start_bias_zvalue_rev, start_bias_pvalue_rev = stats.ranksums(ref_starts_rev, alt_starts_rev)


            if ( start_bias_pvalue_fwd < 0.05 or start_bias_pvalue_fwd < 0.05):
                print str(start_bias_pvalue_fwd) + " " + str(start_bias_pvalue_rev)
                print "%s (%f) vs %s (%f) fails, should it drop %s?" % ( ref_base, ref_freq, alt_base, alt_freq, alt_base)                
#                del base_freqs[ alt_base]
#                exit()

            # Use fishers test to see if there is a strand bias for the SNP
            oddsratio, strand_bias_pvalue = stats.fisher_exact([[len(base_counts[ ref_base ][ 0 ].keys()), 
                                                                 len(base_counts[ ref_base ][ 1 ].keys())], 
                                                                [len(base_counts[ alt_base ][ 0 ].keys()), 
                                                                 len(base_counts[ alt_base ][ 1 ].keys())]])
            if ( strand_bias_pvalue < 0.05 ):
#                pp.pprint( base_counts)
                print "\t".join([str(len(base_counts[ ref_base ][ 0 ].keys())),
                                 str(len(base_counts[ ref_base ][ 1 ].keys())),
                                 str(len(base_counts[ alt_base ][ 0 ].keys())),
                                 str(len(base_counts[ alt_base ][ 1 ].keys()))])

#                print "%d -> %d [%s %s] [%f %f]" % (0, i +1 , ref_base, alt_base, ref_freq, alt_freq)
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
    added_D = 0
    added_I = 0
    for AA in AAs_by_freq:
        AA_line.append( "%s:%.2f%%" % (AA, AA_freqs[ AA ]))

        if ( re.match('insertion', AA ) and not added_D ):
            added_D = 1
            Stanford_format[gene_name].append("%s%d%s" % (ref_AA, AA_number, "i"))
        elif ( re.match('deletion', AA ) and not added_I ):
            added_I = 1
            Stanford_format[gene_name].append("%s%d%s" % (ref_AA, AA_number, "d"))
        else:
            Stanford_format[gene_name].append("%s%d%s" % (ref_AA, AA_number, AA))

        
    AA_changes.append("\t".join( AA_line ))
#    pp.pprint( AA_counts )
#    exit()
    


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


def readin_regions( regions_file ):

    regions = []

    regions_fh = open( regions_file, 'r')
    for l in regions_fh:
        l = l.strip('\n')
        if (re.search('#', l)):
            continue

        fields = l.split("\t") + [None]*99
        for n in range(2, 5):
            if ( fields[ n ]):
                fields[ n ] = int(fields[ n ])

        regions.append(fields[0:6])

    return regions

def print_tab_output( tab_data_lines = None):
#    return

    global tab_lines

    print "\t".join(["pile.pos",
                 "ref/major base",
                 "alt base",
                 "major count",
                 "ref count",

                 "ref freq",
                 "alt freq",
                 "ref starts",
                 "ref starts fwd",
                 "re starts rev",
                 "alt starts",
                 "alt starts fwd",
                 "alt starts rev",
                         
                 "ref read1",
                 "ref read2",
                 
                 "alt read1",
                 "alt read2"])


    print "\n".join( tab_lines )



def print_fasta_consensus():

            
#### get input filename without the path


    consensus_seq = "".join( fasta_consensus )


    if ( minus_strand ):
#        pass
        consensus_seq = revDNA( consensus_seq )


    consensus_seq = re.sub(r'^N+',"", consensus_seq)
    consensus_seq = re.sub(r'N+$',"", consensus_seq)

    if ( len( consensus_seq ) == 0):
        return 

#    if ( len(consensus_seq) == 0):
#        exit()

    header = ">%s_Consensus" % sample_name

    if ( gene_name ):
        header += " %s" % gene_name


    if ( ORF_start > 0 and ORF_end > -1):
        header += " %d-%d" %(ORF_start, ORF_end)

    if ( fasta_fh ):
        fasta_fh.write( header + "\n" )
        fasta_fh.write( consensus_seq + "\n" )
    else:
        print header
        print consensus_seq



def print_codon_changes():


    if ( len( AA_changes) == 0):
        return

    block = gene_name +"\n"
    block +=  "\n".join( AA_changes )
    block += "\n\n"
    block += " ".join(Stanford_format[ gene_name]) + "\n"

    if ( codon_fh ):
        codon_fh.write( block + "\n")
    else:
        print block
    

def readin_bamfile( chrom = "", start = -1, end = -1 ):

    deletion_skipping = 0

    max_depth = 0
    codons    = dict()
    for pile in bam.pileup(chrom, start, end, max_depth=DOWNSAMPLE_DEPTH):

        if ( start > -1 and end > -1) :
            if (pile.pos < start or pile.pos > end):
                continue

        passed_bases = 0

        base_counts  = dict()
        AA_counts = dict()
        insertion_count = 0
        deletion_count  = 0
        in_coding_frame = 0
        AA_number       = 0
        ref_AA          = ""
        ref_base        = 'N'
        if ( fasta_ref ):
            ref_base        = fasta_ref.fetch(str(bam.getrname(pile.tid)), pile.pos, pile.pos + 1 )


        if (max_depth < pile.n ):
            max_depth = pile.n


    # We are in a coding frame...
        if ( output_format & CODON_OUT and (ORF_start - 1  - pile.pos )%3 == 0):
            
            in_coding_frame = 1
            
            (AA_counts, codons) = codons2AA_counts( codons )
            AA_number    = (1+ (pile.pos + 1 - ORF_start - 3)/3)
            genome_pos   = pile.pos+1 - 3;

            if ( first_codon_to_report > 0 and first_codon_to_report > AA_number):
#            print "!!!!!!!" + str(first_codon_to_report) + "\t" + str(codon_number)
                continue

            if ( last_codon_to_report and last_codon_to_report < AA_number):
                continue
            
            if ( fasta_ref ):
                ref_codon = fasta_ref.fetch(str(bam.getrname(pile.tid)), pile.pos - 3, pile.pos )
                ref_AA    = codon2AA( ref_codon )

            find_significant_AAs( AA_counts, genome_pos, ref_AA, AA_number )


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
#                alt = fasta_ref.fetch(str(bam.getrname(pile.tid)), pile.pos, pile.pos+abs(read.indel)+1 )
                alt = read.alignment.seq[ read.qpos ] + "-" * abs(read.indel)
                avg_base_qual = ord( read.alignment.qual[ read.qpos] ) - 33

            else:
                avg_base_qual = ord( read.alignment.qual[ read.qpos] ) - 33
                alt = read.alignment.seq[ read.qpos];

            if ( avg_base_qual > MIN_BASEQ):
                base_counts = update_base_counts(base_counts, alt, 
                                                 avg_base_qual, read.alignment.is_reverse, read.alignment.aend, read_nr)

                if ( output_format & CODON_OUT ):
                    codons = grow_codons( codons, read.alignment.qname, alt, avg_base_qual, 
                                          read.alignment.is_reverse, read.alignment.aend, read_nr, in_coding_frame)



        if ( deletion_skipping > 0 ):
            deletion_skipping -= 1
        else:
            consensus_base = find_consensus_base( base_counts, pile.pos, ref_base )
            # If the consensus base contains a deletion, set the level of skipping here.
            deletion_skipping =  consensus_base.count('-')


            fasta_consensus.append( consensus_base )


    if ( output_format & TAB_OUT ):   
        print_tab_output()

    if ( output_format & FASTA_OUT ):   
        print_fasta_consensus()


    if ( output_format & CODON_OUT ):   
        print_codon_changes()


def get_and_parse_options():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-b', '--bam')
    parser.add_argument('-R', '--reference')

    parser.add_argument('-f', '--fasta_output', action="store_true")
    parser.add_argument('-c', '--codon_output', action="store_true")
    parser.add_argument('-B', '--both_output',  action="store_true")
    parser.add_argument('-T', '--tab_output',   action="store_true" )

    parser.add_argument('-p', '--prefix_output', dest='prefix')
    parser.add_argument('-r', '--regions_file', dest='regions')
    parser.add_argument('--region')

    args = parser.parse_args()

    if (args.regions):
        regions = readin_regions( args.regions )

    if (args.region):
        (chromo, start, end) = re.split(':|-', args.region)

        regions = [[None, chromo, int(start), int(end), None, None]]

    output = 0 

    # multiple outputs need a file prefix so we can create multiple files
    if ( args.both_output and not args.prefix ):
        print "Needs a prefix when doing fasta and codon change output"
        exit( -1 )

        
    if ( args.codon_output or args.both_output or args.tab_output):
        if ( not args.reference ):
            print "A reference (-R) is needed when doing codon or tab output"
            exit( -1 )

    if ( args.fasta_output or args.both_output):
        output += FASTA_OUT
        if ( args.prefix ):
            global fasta_fh
            fasta_fh = open( args.prefix + "_consenus.fasta", 'w')


    if ( args.codon_output or args.both_output):
        output += CODON_OUT
        if ( args.prefix ):
            global codon_fh
            codon_fh = open( args.prefix + "_codons.xls", 'w')

    if ( args.tab_output):
        output = TAB_OUT


    if ( not output ):
        print "Select one of the output formats (-f, -c, -B or -T)"
        exit( -1 )

    return (args, regions, output)


# ---------------- MAIN LOOP ------------------------


fasta_fh = None
codon_fh = None

(args, regions, output_format) = get_and_parse_options()


bamfile   = args.bam
bam       = pysam.Samfile( bamfile, "rb" )

fasta_ref = None

if ( args.reference ):
    fasta_ref = args.reference
    fasta_ref = pysam.Fastafile( fasta_ref );



bamfile = bamfile.replace('.bam', '')
b = re.search(".*\/*(.*)$", bamfile)
if ( b.group(1) ):
    sample_name = b.group(1)
else:
    sample_name = bamfile

start        = -1
end          = -1
minus_strand = 0

Stanford_format = dict()
fasta_consensus = []
AA_changes      = []
tab_lines       = []


ORF_start = None
ORF_end   = None
gene_name = None

first_codon_to_report = None
last_codon_to_report  = None

# The user asked us to look at specific regions
if ( len( regions ) >= 1 ):

    for region in regions:
        Stanford_format = dict()
        fasta_consensus = []
        AA_changes      = []
        tab_lines       = []

        minus_strand = 0

        (gene_name, reference, region_start, region_end, first_codon_to_report, last_codon_to_report) = region
        if ( region_start > region_end ):
            minus_strand = 1
            (region_start, region_end) = (region_end, region_start)

        if ( output_format & CODON_OUT or output_format & FASTA_OUT ):
            ORF_start = region_start
            ORF_end   = region_end

        if ( gene_name ):
            Stanford_format[ gene_name ] = []


        readin_bamfile(reference, int( region_start), int( region_end))

# Generic full out consensus calling.
else:
    readin_bamfile()

bam.close()

