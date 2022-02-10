#!/usr/bin/env python

"""
Author:		Pontus Skoglund
Contact: 	pontus.skoglund@gmail.com
Date: 		January 23, 2014
Citation:	P Skoglund, BH Northoff, MV Shunkov, AP Derevianko, S Paabo, J Krause, M Jakobsson (2014) Separating endogenous ancient DNA from modern day contamination in a Siberian Neandertal, PNAS, advance online 27 January

Usage:		python pmdtools.py <SAM formatted data with MD field present from stdin> [options]

Example:	#remove all sequence reads with PMD score <3
		samtools view -q 30 mybamfile.bam | python pmdtools.py --header --threshold 3 > samtools view -Sb - > mybamfile_filtered.bam

		#for more options:
		python pmdtools.py --help

		(for specification on the SAM format and a the samtools suite, see Li, Handsaker et al. 2009, Bioinformatics)
"""

import sys
from optparse import OptionParser
import math
import subprocess


usage = "usage: python %prog <SAM formatted data with MD field present from stdin> [options] "
parser = OptionParser(usage=usage, version="%prog v0.50")


parser.add_option("--adjustss", action="store_true", dest="adjustss",help="strand orientation aware penalization of base qualities for single stranded libraries (use with --CpG for udg-treated ss libs)",default=False)

parser.add_option("-n", "--number", action="store", type="int", dest="maxreads",help="stop after these many reads have been processed",default=(10**20))
parser.add_option("-c", "--chromosome", action="store", type="string", dest="chromosome",help="only process data from this chromosome",default=False)
parser.add_option("-m", "--requiremapq", action="store", type="int", dest="mapq",help="only process sequences with mapping quality at least this great",default=0)
parser.add_option("--readlength", action="store", type="int", dest="readlength",help="only process sequences with this read length",default=0)
parser.add_option("--maxlength", action="store", type="int", dest="maxlength",help="only process sequences with max this read length",default=300)
parser.add_option("--minlength", action="store", type="int", dest="minlength",help="only process sequences with min this read length",default=0)
parser.add_option("--maxGC", action="store", type="float", dest="maxGC",help="only process sequences with max this GC content of the aligning reference sequence",default=1.0)
parser.add_option("--minGC", action="store", type="float", dest="minGC",help="only process sequences with min this GC content of the aligning reference sequence",default=0.0)
parser.add_option("-q", "--requirebaseq", action="store", type="int", dest="baseq",help="only process bases with base quality at least this great",default=0)
parser.add_option("-d", "--deamination", action="store_true", dest="deamination",help="output base frequencies in the read at positions where there are C or G in the reference",default=False)
parser.add_option("--CpG","--UDGplus", action="store_true", dest="cpg",help="only use Cs and Gs in CpG context",default=False)

parser.add_option("--ss", action="store_true", dest="ss",help="single stranded",default=False)

parser.add_option("--PMDSprim", action="store_true", dest="PMDSprim",help="PMDSprim",default=False)
parser.add_option("--PMDSprimthreshold", action="store",type="float", dest="PMDSprimthreshold",help="PMDSprimthreshold",default=False)
parser.add_option("--UDGhalf", action="store_true", dest="UDGhalf",help="only use Cs and Gs in CpG context, the first and last base are used regardless of dinucleotide context",default=False)
#parser.add_option("--UDGplus", action="store_true", dest="UDGplus",help="only use Cs and Gs in CpG context (synonymous with option --CpG)",default=False)
parser.add_option("--UDGminus", action="store_true", dest="UDGminus",help="use all bases (placeholder)",default=False)


parser.add_option("--EcoliCpG", action="store_true", dest="EcoliCpG",help="Ecoli CpG (first 5' position and two 3' positions and CpG context)",default=False)
parser.add_option("--Ecoli", action="store_true", dest="Ecoli",help="Ecoli (first 5' position and two 3' positions only)",default=False)


parser.add_option("--maskss", action="store_true", dest="maskss",help="mask to lower case bases where ref is C (taking strand into account)",default=False)
parser.add_option("--Leipzigsimple", action="store_true", dest="Leipzigsimple",help="Leipzigsimple",default=False)


parser.add_option("--customterminus", action="store",type="string", dest="customterminus",help="customterminus",default=False)

parser.add_option("--flagss", action="store_true", dest="flagss",help="flag SNPs (taking strand into account)",default=False)

parser.add_option("--noCpG", action="store_true", dest="nocpg",help="dont use Cs and Gs in CpG context",default=False)
parser.add_option("--first", action="store_true", dest="first",help="outputs the deamination rate at the first position only, but with a standard error",default=False)
parser.add_option("--range", action="store", type="int", dest="range",help="output deamination patterns for this many positions from the sequence terminus (default=30)",default=30)
parser.add_option("--polymorphism_ancient", action="store", type="float", dest="polymorphism_ancient",help="True biological polymorphism between the ancient individual and the reference sequence",default=0.001)
parser.add_option("--polymorphism_contamination", action="store", type="float", dest="polymorphism_contamination",help="True biological polymorphism between the contaminants and the reference sequence",default=0.001)
parser.add_option("--PMDpparam", action="store", type="float", dest="PMDpparam",help="parameter p in geometric probability distribution of PMD",default=0.3)
parser.add_option("--PMDconstant", action="store", type="float", dest="PMDconstant",help="constant C in geometric probability distribution of PMD",default=0.01)
parser.add_option("--noclips", action="store_true", dest="noclips",help="no clips",default=False)
parser.add_option("--noindels", action="store_true", dest="noindels",help="no indels",default=False)
parser.add_option("--onlyclips", action="store_true", dest="onlyclips",help="only clips",default=False)
parser.add_option("--onlydeletions", action="store_true", dest="onlydeletions",help="only deletions",default=False)
parser.add_option("--onlyinsertions", action="store_true", dest="onlyinsertions",help="only insertions",default=False)
parser.add_option("--nodeletions", action="store_true", dest="nodeletions",help="no deletions",default=False)
parser.add_option("--noinsertions", action="store_true", dest="noinsertions",help="no insertions",default=False)
parser.add_option("--notreverse", action="store_true", dest="notreverse",help="no reverse complement alignments",default=False)
parser.add_option("-p", "--printDS", action="store_true", dest="printDS",help="print PMD scores",default=False)
parser.add_option("--verbose", action="store_true", dest="verbose",help="verbose",default=False)
parser.add_option("--printalignments", action="store_true", dest="printalignments",help="print human readable alignments",default=False)

parser.add_option("--maskterminaldeaminations", action="store",type="int", dest="maskterminaldeaminations",help="mask terminal deaminations",default=False)
parser.add_option("--maskterminalbases", action="store",type="int", dest="maskterminalbases",help="mask terminal bases",default=False)


parser.add_option("-t", "--threshold", type="float", dest="threshold",help="only output sequences with PMD score above this threshold",default=(-20000.0))
parser.add_option("--upperthreshold", type="float", dest="upperthreshold",help="only output sequences with PMD score below this threshold",default=(1000000.0))
parser.add_option("--perc_identity", type="float", dest="perc_identity",help="only output sequences with percent identity above this threshold",default=0.0)
parser.add_option("-a", "--adjustbaseq", action="store_true", dest="adjustbaseq",help="apply PMD-aware adjustment of base quality scores specific to C>T and G>A mismatches to the reference",default=False)
parser.add_option("--adjustbaseq_all", action="store_true", dest="adjustbaseq_all",help="apply PMD-aware adjustment of base quality scores regardless of observed bases",default=False)
parser.add_option("--dry", action="store_true", dest="dry",help="print SAM input without any filters",default=False)
parser.add_option("--samtoolspath", action="store", dest="samtoolspath",help="full path to samtools",default='samtools')
parser.add_option("--basecomposition", action="store_true", dest="basecomposition",help="basecomposition",default=False)
parser.add_option("-r","--refseq", action="store", dest="refseq",help="refseq",default=False)
parser.add_option("--header", action="store_true", dest="header",help="output the SAM header",default=False)
parser.add_option("--estimate", action="store_true", dest="estimate",help="two-terminus estimate of contamination",default=False)
parser.add_option("--estimatebase", action="store",type="int", dest="estimatebase",help="position of base used fortwo-terminus estimate of contamination",default=0)
parser.add_option("--platypus", action="store_true", dest="platypus",help="output big list of base frequencies for platypus",default=False)
parser.add_option("--writesamfield", action="store_true", dest="writesamfield",help="add 'DS:Z:<PMDS>' field to SAM output, will overwrite if already present",default=False)
parser.add_option("-b", "--basic", action="store", type="int", dest="basic",help="only output reads with a C>T mismatch this many basepairs from the 5' end",default=0)
parser.add_option("--terminal", action="store_true", dest="terminal",help="only output reads with a C>T mismatch in either 3' or 5' end",default=False)
parser.add_option("--stats", action="store_true", dest="stats",help="output summarizing statistics to stderr",default=False)
(options, args) = parser.parse_args()

#if options.UDGplus:
#	options.CpG=True

def translate(inbase):
	if inbase == 'A': outbase = 'T'
	elif inbase == 'T': outbase = 'A'
	elif inbase == 'G': outbase = 'C'
	elif inbase == 'C': outbase = 'G'
	elif inbase == 'N': outbase = 'N'	
	elif inbase == 'Y': outbase = 'Y'
	elif inbase == 'R': outbase = 'R'
	elif inbase == '-': newseq += '-'
	return outbase

def revcomp(inseq):
	inseq= inseq[::-1]#.upper()
	newseq=''
	for inbase in inseq:
		if inbase == 'A': newseq += 'T'
		elif inbase == 'T': newseq += 'A'
		elif inbase == 'G': newseq += 'C'
		elif inbase == 'C': newseq += 'G'
		elif inbase == 'Y': newseq += 'Y'
		elif inbase == 'R': newseq += 'R'
		elif inbase == 'a': newseq += 't'
		elif inbase == 't': newseq += 'a'
		elif inbase == 'g': newseq += 'c'
		elif inbase == 'c': newseq += 'g'
		elif inbase == 'N': newseq += 'N'	
		elif inbase == '-': newseq += '-'
	return newseq

def phred2prob(Q):
	return 10.0 ** (-Q/10.0)

def prob2phred(P):
	return -10.0*math.log(P,10)


def L_match(fposition,fmodel,fquals,fpoly):
	P_damage= 	float(fmodel[fposition]) 
	P_error= 	phred2prob( (ord(fquals[fposition])-33))/3.0
	P_poly=		fpoly
	P_match= 	(1.0-P_damage)*(1.0-P_error)*(1.0-P_poly) + (P_damage*P_error*(1.0-P_poly)) + (P_error*P_poly * (1.0-P_damage))

	return P_match


def L_mismatch(fposition,fmodel,fquals,fpoly):
	P_damage= 	float(fmodel[fposition]) 
	P_error= 	phred2prob( (ord(fquals[fposition])-33))/3.0
	P_poly=		fpoly
	P_match= 	(1.0-P_damage)*(1.0-P_error)*(1.0-P_poly) + (P_damage*P_error*(1.0-P_poly)) + (P_error*P_poly * (1.0-P_damage))
	P_mismatch=	1.0-P_match
	return P_mismatch
	
	

	
def L_match_SS(fposition,zposition,fmodel,fquals,fpoly):
	P_damage= 	float(fmodel[fposition]) 
	P_damage2= 	float(fmodel[zposition]) 
	P_error= 	phred2prob( (ord(fquals[fposition])-33))/3.0
	P_poly=		fpoly
	P_match= 	(1.0-P_damage)*(1.0-P_error)*(1.0-P_poly)*(1.0-P_damage2)  + (P_damage*P_error*(1.0-P_poly)*(1.0-P_damage2) ) + (P_damage2*P_error*(1.0-P_poly)*(1.0-P_damage) ) + (P_error*P_poly*(1.0-P_damage2) * (1.0-P_damage))

	return P_match

def L_mismatch_SS(fposition,zposition,fmodel,fquals,fpoly):
	P_damage= 	float(fmodel[fposition]) 
	P_damage2= 	float(fmodel[zposition]) 
	P_error= 	phred2prob( (ord(fquals[fposition])-33))/3.0
	P_poly=		fpoly
	P_match= 	(1.0-P_damage)*(1.0-P_error)*(1.0-P_poly)*(1.0-P_damage2)  + (P_damage*P_error*(1.0-P_poly)*(1.0-P_damage2) ) + (P_damage2*P_error*(1.0-P_poly)*(1.0-P_damage) ) + (P_error*P_poly*(1.0-P_damage2) * (1.0-P_damage))
	P_mismatch=	1.0-P_match
	return P_mismatch



def Newbaseq(fposition,fmodel,fquals):
	P_damage= 	float(fmodel[fposition]) 
	P_error= 	phred2prob( (ord(fquals[fposition])-33))/3.0
	NewErrorP= 1.0 - ((1.0-P_damage) * (1.0-P_error))
	return NewErrorP

def geometric(pval,kval,constant):
	return ((1.0-pval)**(kval-1))*pval + constant

FNULL = open('/dev/null', 'w')
def fa_get (ffile,chrom,fstart,fend):
	chrom=str(chrom)#'chr'+str(chrom)
	location=str(chrom)+':'+str(fstart)+'-'+str(fend)
	cmd_line=[options.samtoolspath, 'faidx', ffile,location]
	outp_file=subprocess.Popen(cmd_line, stdout=subprocess.PIPE,stderr=FNULL)
	pileupline = outp_file.stdout.read().split()
	pileupline =''.join(pileupline[1:])
	if len(pileupline)<1:
			print 'no such reference sequence',cmd_line
	if options.verbose:
		print ' '.join(cmd_line)
	return pileupline
	
		

maxlen=options.maxlength
 
modern_model_deam='0.001'
modern_model_deam=modern_model_deam.split('\n')*1000 
modern_model_deam=[float(l) for l in modern_model_deam]

ancient_model_deam=[geometric(options.PMDpparam,l,options.PMDconstant) for l in range(1,1000)]

if options.adjustbaseq_all:
	adjustment_model_deam=[geometric(options.PMDpparam,l,0.0) for l in range(1,1000)] ###constant is 0.0 here in contrast to the model used to compute PMD scores

start_dict= {}

###base composition
start_count = 0
rev_start_count = 0
not_counted = 0
imperfect = 0

mismatch_dict={}


import re
cigarparser = re.compile("([0-9]*)([A-Z])")


start_dict_rev= {}
mismatch_dict_rev={}


mismatch_dict_CpG={}
mismatch_dict_CpG_rev={}
mismatch_dict_CpG_rev_ss={}

firstC=0
firstT=0


clipexcluded=0
indelexcluded=0
noMD=0
noGCexcluded=0
excluded_threshold=0
passed=0
noquals=0
maskings=0


CCandCC=0
CTandCC=0
CCandCT=0
CTandCT=0


estimator_list=[]


composition_dict={}
composition_dict_rev={}

line_counter = 0
for line in sys.stdin: 
	if '@' in line[0]: 
		if options.header:
			print line.rstrip('\n')
		continue
	line_counter +=1
	line=line.rstrip('\n')
	col = line.split('\t')
	readname = col[0]
	position = int(col[3])
	chromosome = col[2]

	if options.chromosome:
		if chromosome != options.chromosome:continue
	MAPQ = int(col[4])
	read = col[9]
	readlen = len(read)
	quals= col[10]
	flag = col[1]
	position = int(col[3])
	cigar=col[5]

	if len(quals) <2:
		noquals+=1
		continue

	if options.noinsertions:
		if 'I' in cigar:continue 
	if options.nodeletions:
		if 'D' in cigar:continue 
	if options.onlyinsertions:
		if 'I' not in cigar:continue 
	if options.onlydeletions:
		if 'D' not in cigar:continue 
	if options.noindels:
		if 'I' in cigar or 'D' in cigar:
			indelexcluded +=1
			continue
	if options.noclips:
		if 'S' in cigar or 'H' in cigar or 'N' in cigar or 'P' in cigar:
			clipexcluded +=1
			continue
	if options.onlyclips:
		if 'S' not in cigar:
			continue
	if 'H' in cigar or 'P' in cigar or 'N' in cigar:
		print >>sys.stderr,'cigar found:',cigar,'PMDtools only supports cigar operations M, I, S and D, the alignment has been excluded'
		continue
	if MAPQ < options.mapq: 
		continue

	if options.readlength >0: 
		if options.readlength != len(read):continue
	if options.minlength >0: 
		if options.minlength > len(read):continue
	if options.maxlength != 300: 
		if options.maxlength < len(read):continue
	if options.chromosome:
		if chromosome != options.chromosome: continue



	if flag.isdigit():
		if int(flag) & 16:
			reverse = True
		else:
			reverse = False
	else:
		if 'r' in flag: reverse = True
		elif 'r' not in flag: reverse = False

	if options.notreverse:
		if reverse: continue

	DSfield=False
	if 'DS:Z:' in line:
		DSfield=True
		PMDS= float(line.split('DS:Z:')[1].rstrip('\n').split()[0])
		LR=PMDS
		#print PMDS
	
	


	"""
	Recreate reference sequence from MD field
	"""
	if (DSfield == False) or (options.writesamfield) or (options.basic > 0) or (options.terminal) or (options.perc_identity > 0.01) or (options.printalignments) or options.adjustss or (options.adjustbaseq) or (options.adjustbaseq_all) or options.deamination or options.dry  or options.estimate or options.first:
		
		read=col[9]
		ref_seq=''
		newread=''
	
		import re
		try:
			MD=line.split('MD:Z:')[1].split()[0].rstrip('\n')
		except IndexError:
			noMD+=1
			continue

		MDlist=re.findall('(\d+|\D+)',MD)


		MDcounter=0
		alignment=''
		for e in MDlist:
			if e.isdigit():
				e=int(e)
				alignment += '.'*e
				ref_seq+= read[MDcounter:MDcounter+e]
				newread+= read[MDcounter:MDcounter+e]
				MDcounter+= int(e)

			elif '^' in e:
				ef=e.lstrip('^')
				alignment += ef
				continue
			elif e.isalpha():
				alignment += e
				ref_seq+= e
				newread+= read[MDcounter]
				MDcounter+=len(e)

		if 'I' in cigar or 'S' in cigar:

			# find insertions and clips in cigar
			insertions=[]
			deletions=[]
			softclips=[]
			hardclips=[]
			paddings=[]
			cigarcount=0
			cigarcomponents=cigarparser.findall(cigar)
			for p in cigarcomponents:
				cigaraddition=int(p[0])
				if 'I' in p[1]: 
					for c in range(cigarcount,cigarcount+cigaraddition): insertions.append(c)
				elif  'S' in p[1]:
					for c in range(cigarcount,cigarcount+cigaraddition): softclips.append(c)
				cigarcount += int(p[0])
			# end cigar parsing

			# redo the read and ref using indel and clip info
			ref_seq=''
			newread =''
			alignmentcounter=0
			for x,r in zip(xrange(0,len(col[9])),read):
				if x in insertions:
					ref_seq += '-'
					newread += read[x]
				elif x in softclips:
					ref_seq += '-'
					newread += read[x]
				else:
					newread += read[x]
					refbasealn=alignment[alignmentcounter]
					if refbasealn == '.':
						ref_seq += read[x]
					else:
						ref_seq += refbasealn
					alignmentcounter +=1
					
		if reverse:
			read = revcomp(read)
			ref_seq = revcomp(ref_seq)
			quals = quals[::-1]
		real_read=read
		real_ref_seq=ref_seq
		
		
	if options.maskterminaldeaminations or options.maskterminalbases:
		maskedseq=real_read
			
	GCcontent=1.0*(real_ref_seq.count('G')+real_ref_seq.count('C'))/readlen
	if GCcontent > options.maxGC:continue
	elif GCcontent < options.minGC:continue
	if 'G' not in real_ref_seq and 'C' not in real_ref_seq and 'T' not in real_ref_seq and 'A' not in real_ref_seq:
		print >>sys.stderr,'bad reference sequence reconstruction:',real_ref_seq
		print >>sys.stderr,'SAM line:',line
		continue
		#exit(1)
	
	if options.basecomposition:
		backoffset=10
		if reverse:
			endpos=position
			startpos=position+len(real_read)
		else:
			startpos=position
			endpos=position+len(real_read)
			
		"""
		5' end
		"""	
		largerefseq=fa_get(options.refseq,chromosome,startpos-backoffset,startpos+options.range)
		largerefseq=largerefseq.upper()
		if len(largerefseq)<1:continue
		#print largerefseq
		#if reverse:
		#	largerefseq=revcomp(largerefseq)
		#print largerefseq,real_read	
		for i in range(-backoffset,options.range):
					#print i+backoffset, len(largerefseq)
					base=largerefseq[min([i+backoffset,len(largerefseq)])]
					thekey='5'+base+str(i)
					if thekey in composition_dict.keys():
						addition = composition_dict[thekey]
						addition += 1
						composition_dict[thekey] = addition
					else:
						composition_dict[thekey] = 1
						
		"""
		3' end
		"""	
		largerefseq=fa_get(options.refseq,chromosome,endpos-options.range,endpos+backoffset)
		largerefseq=largerefseq.upper()[::-1]
		if len(largerefseq)<1:continue
		#if reverse:
		#	largerefseq=revcomp(largerefseq)
			
		for i in range(-backoffset,options.range):
					
					base=largerefseq[min([i+backoffset,len(largerefseq)])]
					thekey='3'+base+str(i)
					if thekey in composition_dict_rev.keys():
						addition = composition_dict_rev[thekey]
						addition += 1
						composition_dict_rev[thekey] = addition
					else:
						composition_dict_rev[thekey] = 1
		#continue
	

	"""
	basic filter
	prints the SAM line if a C>T mismatch with sufficient base quality is observed in the first n bases, where n is specified
	"""
	if options.basic > 0:
		#start_position = len(real_read) - len(real_read.lstrip('-'))
		for a,b,x in zip(real_read,real_ref_seq,range(0,len(real_ref_seq))):

			if a == 'N': break
			elif b == 'N': break
			elif a=='-':break
			elif b=='-':break
			i = x #- start_position
			if i >= readlen: break ###20
			if i >= options.basic: break
			if options.cpg:
				if b == 'C' and a=='T' and ((ord(quals[i])-33) > options.baseq) and (real_ref_seq[i+1] == 'G'): 
					print line.rstrip('\n')
					break

			elif b == 'C' and a=='T' and ((ord(quals[i])-33) > options.baseq): 
				print line.rstrip('\n')
				break
				
	"""
	basic terminal
	prints the SAM line if a C>T mismatch with sufficient base quality is observed in the first and last base
	"""
	if options.terminal:
		a,b = real_read[0],real_ref_seq[0]
		if b == 'C' and a=='T' and ((ord(quals[0])-33) > options.baseq): 
			#print a,b
			print line.rstrip('\n')
			continue
			#break
		a,b = real_read[::-1][0],real_ref_seq[::-1][0]
		#i=len(real_read)-1
		if b == 'G' and a=='A' and ((ord(quals[::-1][0])-33) > options.baseq): 
			print line.rstrip('\n')
			continue
			#break
		#continue

	"""
	first base
	prints the deamination rate at the first base and a standard error computed by jackknife over reads
	"""
	if options.first:
		b=real_ref_seq[0]
		a=real_read[0]
		if options.cpg:
			if b == 'C' and a=='T' and ((ord(quals[0])-33) >= options.baseq) and (real_ref_seq[1] == 'G'): 
				firstT +=1
			if b == 'C' and a=='C' and ((ord(quals[0])-33) >= options.baseq) and (real_ref_seq[1] == 'G'): 
				firstC +=1
		else:
			if b == 'C' and a=='T' and ((ord(quals[0])-33) >= options.baseq): 
				firstT +=1
			if b == 'C' and a=='C' and ((ord(quals[0])-33) >= options.baseq): 
				firstC +=1
				
				
	if options.Leipzigsimple:

		Leipzigsimple=False
		for leipnum in [0,-1,-2]:
			b=real_ref_seq[leipnum]
			a=real_read[leipnum]
			if b == 'C' and a=='T' and ((ord(quals[leipnum])-33) >= options.baseq): 
				Leipzigsimple =True
		if Leipzigsimple:
			line=line.rstrip('\n')+'\t'+'LS:Z:'+'1'
		else:
			line=line.rstrip('\n')+'\t'+'LS:Z:'+'0'
			
	if options.customterminus != False:
		
		Leipzigsimple=False
		customterminus_list=[int(l) for l in options.customterminus.split(',')]
		if len(real_ref_seq) < max(customterminus_list) or len(real_ref_seq) < abs(min(customterminus_list)):
			continue
		
		for leipnum in customterminus_list:
			b=real_ref_seq[leipnum]
			a=real_read[leipnum]
			if leipnum <0:
					if options.ss:
						if b == 'C' and a=='T' and ((ord(quals[leipnum])-33) >= options.baseq): 
							Leipzigsimple =True
					else:
						if b == 'G' and a=='A' and ((ord(quals[leipnum])-33) >= options.baseq): 
							Leipzigsimple =True
			else:
				if b == 'C' and a=='T' and ((ord(quals[leipnum])-33) >= options.baseq): 
					Leipzigsimple =True
		if Leipzigsimple:
			print line
		else:
			line=line.rstrip('\n')+'\t'+'LS:Z:'+'0'
			continue

			
	if options.perc_identity > 0.01 or options.printalignments:
		"""
		divergence filter
		"""
		match=0
		mismatch=0
		mismatch_string=''
		for a,b in zip(real_read,real_ref_seq):
			thesebases=[a,b]
			if '-' in thesebases: 
				mismatch_string+='-'
				continue
			if a == 'N': continue
			if b == 'N': continue
			if a == b: 
				mismatch_string+='|'
				match +=1
			elif a!=b: 
				mismatch_string+='x'
				if 'C' == b and 'T' == a: continue
				if 'G' == b and 'A' == a: continue
				mismatch +=1	
		
		try:
			perc_identity=1.0*match/(match+mismatch)
		except ZeroDivisionError:
			continue
		if perc_identity < options.perc_identity:continue





	"""
	start PMD score computations
	"""

	if (DSfield == False) or (DSfield == True and options.writesamfield == True) or (options.basic > 0) or options.adjustss or options.adjustbaseq or options.adjustbaseq_all or options.deamination or options.dry:
		L_D=1.0
		L_M=1.0
		
		L_D_max=1.0
		L_M_max=1.0
		
		L_D_min=1.0
		L_M_min=1.0
		
		"""
		Check for informative sites (does nothing at the moment)
		"""
		#print real_ref_seq[-1]
		if options.cpg:
			if 'CG' not in real_ref_seq:
				#print 'no informative'
				L_D=1.0
				L_M=1.0
		elif options.UDGhalf:
			if 'CG' not in real_ref_seq and real_ref_seq[0] !='C' and real_ref_seq[-1] !='G':
				#print 'no informative'
				L_D=1.0
				L_M=1.0
		else:
			if 'C' not in real_ref_seq and 'G' not in real_ref_seq:
				#print 'no informative'
				L_D=1.0
				L_M=1.0

		newquals=quals
		start_position = 0
		back_start_position = len(real_read)-1
		for a,b,x in zip(real_read,real_ref_seq,range(0,len(real_ref_seq))):
			if 'N' in [a,b]: continue
			i = x - start_position
			z = back_start_position - x 
			qualsrev=quals[::-1]
			if i >= readlen: break ###20

			if options.adjustbaseq_all:
				newprob= adjustment_model_deam[i]+adjustment_model_deam[z] + phred2prob(ord(quals[i])-33)
				#newprob=min(newprob,1.0)
				newphred=int(prob2phred(newprob))
				newqual=chr(int(newphred)+33)						
				newquals=quals[0:i]+newqual+quals[(i+1):]
				
			if options.adjustss:
				if b == 'C' and reverse==False:
					if options.cpg:
						if i+1 >= readlen: break
						if real_ref_seq[i+1] != 'G': continue
					#quality is set to 2
					newqual=chr(2+33)						
					newquals=quals[0:i]+newqual+quals[(i+1):]
						
				elif b == 'A' and reverse==True:
					if options.cpg:
						if i-1 >= readlen: break
						if real_ref_seq[i-1] != 'C': continue
					#quality is set to 2
					newqual=chr(2+33)						
					newquals=quals[0:i]+newqual+quals[(i+1):]

			if (ord(quals[i])-33) < options.baseq:
				#make sure that quality is adjusted even if baseq is below threshold
				if options.adjustbaseq:
					if b == 'C' and a == 'T':
						if options.cpg:
							if i+1 >= readlen: break
							if real_ref_seq[i+1] != 'G': continue
						elif options.nocpg:
							if i+1 >= readlen: break
							if real_ref_seq[i+1] == 'G': continue
						elif options.UDGhalf:
							if i+1 >= readlen: break
							if real_ref_seq[i+1] != 'G' and i != 0: continue
						newprob= Newbaseq(i,ancient_model_deam,quals)
						newphred=int(prob2phred(newprob))
						newqual=chr(int(newphred)+33)						
						newquals=quals[0:i]+newqual+quals[(i+1):]

					if (b == 'G' and a == 'A') or (options.ss and b == 'C' and a == 'T'):
						if options.cpg:
							if i-1 >= readlen: break
							if real_ref_seq[i-1] != 'C': continue
						elif options.nocpg:
							if i-1 >= readlen: break
							if real_ref_seq[i-1] == 'C': continue
						elif options.UDGhalf:
							if real_ref_seq[i-1] != 'C' and z != 0: continue
						newprob= Newbaseq(z,ancient_model_deam,qualsrev)
						newphred=int(prob2phred(newprob))
						newqual=chr(int(newphred)+33)						
						newquals=quals[0:i]+newqual+quals[(i+1):]
						
				continue
				

			"""
			platypus
			"""
			if options.platypus:
				#count from 5' end
				cpgcheck=False
				if i+1 <= len(real_ref_seq):
					if real_ref_seq[i:i+2] == 'CG':
						cpgcheck=True
				#if (i-1 > -1):
				#	if real_ref_seq[i-1:i+1] == 'CG':
				#		cpgcheck=True
				if cpgcheck==True: 
					thekey=b+a+str(i)
					if thekey in mismatch_dict_CpG.keys():
						addition = mismatch_dict_CpG[thekey]
						addition += 1
						mismatch_dict_CpG[thekey] = addition
					else:
						mismatch_dict_CpG[thekey] = 1
				else:
					thekey=b+a+str(i)
					if thekey in mismatch_dict.keys():
						addition = mismatch_dict[thekey]
						addition += 1
						mismatch_dict[thekey] = addition
					else:
						mismatch_dict[thekey] = 1
				
				#count from 3' end
				cpgcheck=False
				if (i-1 > -1):
					if real_ref_seq[i-1:i+1] == 'CG':
						cpgcheck=True
				#if i+1 <= len(real_ref_seq):
				#	if real_ref_seq[i:i+2] == 'CG':
				#		cpgcheck=True
				if cpgcheck==True: 
					thekey=b+a+str(z)
					if thekey in mismatch_dict_CpG_rev.keys():
						addition = mismatch_dict_CpG_rev[thekey]
						addition += 1
						mismatch_dict_CpG_rev[thekey] = addition
					else:
						mismatch_dict_CpG_rev[thekey] = 1
				else:
					thekey=b+a+str(z)
					if thekey in mismatch_dict_rev.keys():
						addition = mismatch_dict_rev[thekey]
						addition += 1
						mismatch_dict_rev[thekey] = addition
					else:
						mismatch_dict_rev[thekey] = 1

						
						
			if options.deamination:
				if b == 'C':
					if options.cpg:
						if i+1 >= readlen: break
						if real_read[i+1] != 'G': continue
					elif options.nocpg:
						if i+1 >= readlen: break
						if real_read[i+1] == 'G': continue
					elif options.UDGhalf:
						if i+1 >= readlen: break
						if real_ref_seq[i+1] != 'G' and i != 0: continue
						
					thekey=b+a+str(i)
					if thekey in mismatch_dict.keys():
						addition = mismatch_dict[thekey]
						addition += 1
						mismatch_dict[thekey] = addition
					else:
						mismatch_dict[thekey] = 1
					

				if b == 'G':
					if options.cpg:
						if i-1 >= readlen: break
						if real_ref_seq[i-1] != 'C': continue
					elif options.nocpg:
						if i-1 >= readlen: break
						if real_ref_seq[i-1] == 'C': continue
					elif options.UDGhalf:
						if i-1 >= readlen: break
						if real_ref_seq[i-1] != 'C' and z != 0: continue
					thekey=b+a+str(z)
					if thekey in mismatch_dict_rev.keys():
						addition = mismatch_dict_rev[thekey]
						addition += 1
						mismatch_dict_rev[thekey] = addition
					else:
						mismatch_dict_rev[thekey] = 1
				continue


			"""	
			compute degradation score
			"""
			if True:
				if i >= readlen:continue
				if b == 'C':
				
					if options.cpg:
						if i+1 >= readlen: break
						if real_ref_seq[i+1] != 'G': continue
					elif options.cpg:
						if i+1 >= readlen: break
						if real_ref_seq[i+1] == 'G': continue
					elif options.UDGhalf:
						if i+1 >= readlen: break
						if real_ref_seq[i+1] != 'G' and i != 0: continue
						
					elif options.EcoliCpG or ('UDG' in line and 'noUDG' not in line):
						#print 'RGdetected',line
						if i+1 >= readlen: break
						if real_ref_seq[i+1] != 'G' and i !=0 and z not in [0,1]: continue

					elif options.Ecoli:
						if i+1 >= readlen: break
						if i !=0 and z not in [0,1]: continue
						
					if a=='T': 
						L_D = L_D * L_mismatch(i,ancient_model_deam,quals,options.polymorphism_ancient) 
						L_M = L_M * L_mismatch(i,modern_model_deam,quals,options.polymorphism_contamination) 
						
						if options.ss:
							L_D = L_D * L_mismatch_SS(i,z,ancient_model_deam,quals,options.polymorphism_ancient) 
							L_M = L_M * L_mismatch_SS(i,z,modern_model_deam,quals,options.polymorphism_contamination) 
						
						
						
						
						
						
						if options.adjustbaseq:
							newprob= Newbaseq(i,ancient_model_deam,quals)
							newphred=int(prob2phred(newprob))
							newqual=chr(int(newphred)+33)
							"""
							print phred2prob(ord(quals[i])),newprob
							print ord(quals[i]),newphred
							print quals[i],newqual
							print quals[0:i],quals[i],quals[(i+1):]
							print quals	
							"""					
							quals=quals[0:i]+newqual+quals[(i+1):]
							
						if options.maskterminaldeaminations != False:
							if (i <= options.maskterminaldeaminations):
								maskedseq=real_read[0:i]+'N'+real_read[(i+1):]
							elif (z <= options.maskterminaldeaminations) and options.ss:
								maskedseq=real_read[0:i]+'N'+real_read[(i+1):]


					elif a=='C': 
						L_D = L_D * L_match(i,ancient_model_deam,quals,options.polymorphism_ancient) 
						L_M = L_M * L_match(i,modern_model_deam,quals,options.polymorphism_contamination) 
						
						if options.ss:
							L_D = L_D * L_match_SS(i,z,ancient_model_deam,quals,options.polymorphism_ancient) 
							L_M = L_M * L_match_SS(i,z,modern_model_deam,quals,options.polymorphism_contamination) 
						
						
							
					if options.PMDSprim and a in ['C','T']:
						L_D_max=L_D_max * L_mismatch(i,ancient_model_deam,quals,options.polymorphism_ancient) 
						L_M_max=L_M_max * L_mismatch(i,modern_model_deam,quals,options.polymorphism_ancient) 
						
						L_D_min=L_D_min * L_match(i,ancient_model_deam,quals,options.polymorphism_ancient) 
						L_M_min=L_M_min * L_match(i,modern_model_deam,quals,options.polymorphism_ancient) 

							

				if b == 'G' and (options.ss ==False):
					if options.cpg:
						if real_ref_seq[i-1] != 'C': continue
					elif options.nocpg:
						if real_ref_seq[i-1] == 'C': continue
					elif options.UDGhalf:
						if real_ref_seq[i-1] != 'C' and z != 0: continue
					if a=='A': 
						L_D = L_D * L_mismatch(z,ancient_model_deam,qualsrev,options.polymorphism_ancient) 
						L_M = L_M * L_mismatch(z,modern_model_deam,qualsrev,options.polymorphism_contamination) 

						if options.adjustbaseq:
							newprob= Newbaseq(z,ancient_model_deam,qualsrev)
							
							newphred=int(prob2phred(newprob))
							newqual=chr(int(newphred)+33)						
							newquals=quals[0:i]+newqual+quals[(i+1):]
							
							
					if options.maskterminaldeaminations != False and options.ss ==False:
							if (z <= options.maskterminaldeaminations):
								maskedseq=real_read[0:i]+'N'+real_read[(i+1):]

					elif a=='G': 
						#try:
						L_D = L_D * L_match(z,ancient_model_deam,qualsrev,options.polymorphism_ancient)
						L_M = L_M * L_match(z,modern_model_deam,qualsrev,options.polymorphism_contamination) 

					if options.PMDSprim and a in ['G','A']:
						L_D_max=L_D_max * L_mismatch(z,ancient_model_deam,quals,options.polymorphism_ancient) 
						L_M_max=L_M_max * L_mismatch(z,modern_model_deam,quals,options.polymorphism_ancient) 
						
						L_D_min=L_D_min * L_match(z,ancient_model_deam,qualsrev,options.polymorphism_ancient) 
						L_M_min=L_M_min * L_match(z,modern_model_deam,qualsrev,options.polymorphism_ancient) 


		LR= (math.log(L_D/L_M)  )
		if options.PMDSprim:
			maxPMDSval=(math.log(L_D_max/L_M_max)  )
			maxPMDSval=maxPMDSval/readlen
			if LR >0.0:
				LRnumerator=(math.log(L_D_max/L_M_max)  )
			elif LR <0.0:
				LRnumerator=(math.log(L_D_max/L_M_max)  )
			if options.PMDSprimthreshold:
				if maxPMDSval < options.PMDSprimthreshold:
					continue
				
			if options.printDS:
				print LR,maxPMDSval,maxPMDSval,maxPMDSval*readlen,readlen
			#LR=LR/LRnumerator
		quals=newquals	

	if options.adjustbaseq:
		if reverse:
			qualsp=quals[::-1]
		else:
			qualsp=quals
		line='\t'.join(col[0:10])+'\t'+qualsp+'\t'+'\t'.join(col[11:])

	if options.maskterminaldeaminations or options.maskterminalbases:
		readp=maskedseq
		if reverse:
			readp = revcomp(readp)
		real_read=maskedseq
		line='\t'.join(col[0:9])+'\t'+readp+'\t'+'\t'.join(col[10:])

	"""
	add PMDS tag
	"""
	if options.writesamfield == True:
		# remove DS field if present
		if DSfield==True:
			newline=''
			for col in line.split('\t'):
				if 'DS:Z:' in col:
					continue
				else:
					newline += col+'\t'
			line=newline.rstrip('\t')
		
		line=line.rstrip('\n')+'\t'+'DS:Z:'+str(round(LR,3))
		
	if options.flagss:
		# remove DS field if present
		if reverse and 'R' in MD:
			line=line.rstrip('\n')+'\t'+'DD:Z:'+'1'
		elif reverse==False and 'Y' in MD:
			line=line.rstrip('\n')+'\t'+'DD:Z:'+'1'
		else:
			line=line.rstrip('\n')+'\t'+'DD:Z:'+'0'

	if options.maskss:
		maskedread=''
		for a,b,x in zip(real_read,real_ref_seq,range(0,len(real_ref_seq))):
			if b=='C':
				maskedread +=a.lower()
			else:
				maskedread +=a
				
		if reverse:
			maskedread=revcomp(maskedread)
		maskedread=''.join(mr for mr in maskedread if mr !='-')
		print real_read
		print real_ref_seq
		print maskedread
		print ''

	if options.printDS:
		print L_D,'\t',L_M,'\t',L_D/L_M,'\t',LR#   ,'\t',readlen,'\t',perc_identity,'\t',perc_identity*(math.log((L_D/L_M)))

	if options.dry:
		if len(line) <1:continue
		print line.rstrip('\n')
		continue
	
	if options.threshold > (-10000) or options.upperthreshold < (1000000): 

		if LR >= options.threshold and LR < options.upperthreshold:
			print line.rstrip('\n')
		else:
			excluded_threshold +=1



	if options.printalignments:
		if options.threshold > (-10000) or options.upperthreshold < (1000000):
			try:
				LR= (math.log(L_D/L_M)  )
			except: continue
			if LR < options.threshold or LR >options.upperthreshold < (1000000):
				continue

		quals1=''
		quals2=''			
		for q in quals:
			qnum=ord(q)-33
			if qnum <10:
				quals1 +='0'
				quals2+=str(qnum)
			else:
				quals1+=str(qnum)[0]
				quals2+=str(qnum)[1]
		#print MD,cigar,reverse
		#print col[9]
		print real_read
		print mismatch_string
		print real_ref_seq
		print quals
		print quals1
		print quals2
		#print col[10]
		print ''
	passed+=1
	if passed >= options.maxreads:break

if options.first:
	n=firstC+firstT
	freq=1.0*firstT/n
	SE=math.sqrt((freq*(1.0-freq))/n)
	if freq==0.0:
		SE='NA'
	print 'C>T_at_1st_position_and_SE:','\t',freq,'\t',SE#,'\t',n,firstC,firstT


if options.stats:
	print >>sys.stderr,'""""""""""""""""""""""""""""""""'
	print >>sys.stderr,'" excluded due to clipping:',clipexcluded
	print >>sys.stderr,'" excluded due to indels:',indelexcluded
	print >>sys.stderr,'" no MD field:',noMD
	print >>sys.stderr,'" no G or C in ref:',noGCexcluded
	print >>sys.stderr,'" total seqs:',passed
	print >>sys.stderr,'" excluded due to PMD score <',str(int(options.threshold))+':',excluded_threshold
	print >>sys.stderr,'" passed seqs:',(passed-excluded_threshold)
	print >>sys.stderr,'""""""""""""""""""""""""""""""""'


if options.deamination:
	if True:
		pairs=['CT','CA','CG','CC','GA','GT','GC','GG']
		itotaldict={}
		ztotaldict={}
		for i in range(0,options.range):
			itotal=0
			ztotal=0
			for p in pairs:
				thekey=p+str(i)
				try:
					itotal += mismatch_dict[thekey]
				except KeyError: pass
				try:
					ztotal += mismatch_dict_rev[thekey]
				except KeyError: pass
			itotaldict[i]=itotal
			ztotaldict[i]=ztotal

		print 'z\t','\t'.join(pairs)


		for i in range(0,options.range):
			print str(i)+'\t',
			for p in pairs:
				thekey=p+str(i)
				if 'C' in p[0]:
					try:
						thecount=mismatch_dict[thekey]
					except KeyError: 
						print '0.00000\t',
						continue
					thetotal=itotaldict[i]
					frac=1.0*thecount/thetotal
				if 'G' in p[0]:
					try:
						thecount=mismatch_dict_rev[thekey]
					except KeyError: 
						print '0.00000\t',
						continue
					thetotal=ztotaldict[i]
					frac=1.0*thecount/thetotal
				print str(round(frac,5))+'\t',
			print ''

if options.platypus:
	if True:
		pairs=['CT','CA','CG','CC','GA','GT','GC','GG','AA','AT','AC','AG','TA','TT','TC','TG'] #
		itotaldict={}
		ztotaldict={}
		CpG_itotaldict={}
		CpG_ztotaldict={}
		
		
		
		for base in ['C','T','G','A']:
			for i in range(0,options.range):
				itotal=0
				ztotal=0
				for base2 in ['C','T','G','A']:
					thekey=base+base2+str(i)
					#print thekey
					try:
						itotal += mismatch_dict[thekey]
					except KeyError: pass
					try:
						ztotal += mismatch_dict_rev[thekey]
					except KeyError: pass
				itotaldict[str(i)+base]=itotal
				ztotaldict[str(i)+base]=ztotal
				#print i,base,itotal

				
		for base in ['C','G','T','A']:
			for i in range(0,options.range):
				itotal=0
				ztotal=0
				for base2 in ['C','T','G','A']:
					thekey=base+base2+str(i)
					#print thekey
					try:
						itotal += mismatch_dict_CpG[thekey]
					except KeyError: pass
					try:
						ztotal += mismatch_dict_CpG_rev[thekey]
					except KeyError: pass
				CpG_itotaldict[str(i)+base]=itotal
				CpG_ztotaldict[str(i)+base]=ztotal
				#print i,base,itotal
				
		print 'z\t',
		for p in pairs:
			print p+'5'+'\t',

		for p in pairs:
			print p+'3'+'\t',
		for p in pairs:
			print p+'_CpG_5'+'\t',

		for p in pairs:
			print p+'_CpG_3'+'\t',
		print ''
		
		
		for i in range(0,options.range):
			print str(i)+'\t',
			for p in pairs:
				thekey=p+str(i)
				try:
					thecount=mismatch_dict[thekey]
				except KeyError: 
					print '0\t',
					continue
				thetotal=itotaldict[str(i)+p[0]]
				frac=1.0*thecount/thetotal
				print str(round(frac,4))+'\t',
				#print str(frac)+'\t',
			#print ''
			
			
			for p in pairs:
				thekey=p+str(i)
				try:
					thecount=mismatch_dict_rev[thekey]
				except KeyError: 
					print '0\t',
					continue
				thetotal=ztotaldict[str(i)+p[0]]
				frac=1.0*thecount/thetotal
				print str(round(frac,4))+'\t',
				#print str(frac)+'\t',

			pairs=['CT','CA','CG','CC','GA','GT','GC','GG','AA','AT','AC','AG','TA','TT','TC','TG']
			for p in pairs:
				thekey=p+str(i)
				try:
					thecount=mismatch_dict_CpG[thekey]
				except KeyError: 
					print '0\t',
					continue
				thetotal=CpG_itotaldict[str(i)+p[0]]
				frac=1.0*thecount/thetotal
				print str(round(frac,4))+'\t',
				#print str(frac)+'\t',
			#print ''
			
			
			for p in pairs:
				thekey=p+str(i)
				try:
					thecount=mismatch_dict_CpG_rev[thekey]
				except KeyError: 
					print '0\t',
					continue
				thetotal=CpG_ztotaldict[str(i)+p[0]]
				frac=1.0*thecount/thetotal
				print str(round(frac,4))+'\t',
				#print str(frac)+'\t',
			print ''
			
			
if options.basecomposition:
	#print composition_dict
	#print composition_dict_rev
	if True:
		pairs=['5T','5A','5G','5C','3T','3A','3G','3C']
		itotaldict={}
		ztotaldict={}
		for i in range(-backoffset,options.range):
			itotal=0
			ztotal=0
			for p in pairs:
				thekey=p+str(i)
				try:
					itotal += composition_dict[thekey]
				except KeyError: pass
				try:
					ztotal += composition_dict_rev[thekey]
				except KeyError: pass
			itotaldict[i]=itotal
			ztotaldict[i]=ztotal

		print 'z\t','\t'.join(pairs)


		for i in range(-backoffset,options.range):
			print str(i)+'\t',
			for p in pairs:
				thekey=p+str(i)
				if '5' in p[0]:
					try:
						thecount=composition_dict[thekey]
					except KeyError: 
						print '0.00000\t',
						continue
					thetotal=itotaldict[i]
					frac=1.0*thecount/thetotal
				if '3' in p[0]:
					try:
						thecount=composition_dict_rev[thekey]
					except KeyError: 
						print '0.00000\t',
						continue
					thetotal=ztotaldict[i]
					frac=1.0*thecount/thetotal
				print str(round(frac,5))+'\t',
			print ''