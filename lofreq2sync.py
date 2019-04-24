#!/usr/bin/python
# ***************************************************************
# Name:      lofreq2sync.py
# Purpose:   From https://code.google.com/archive/p/popoolation2/wikis/Manual.wiki
#	     The program converts lofreq VCF files to the 'synchronized' format. This is 
#            the standard format used with all PoPoolation2 scripts (https://code.google.com/archive/p/popoolation2/wikis/Manual.wiki). 
#	     The synchronized file contain the allele frequencies for all bases in the 
#	     reference genome and for all populations being analyzed. From the above link, following an example of two populations:
#
#	     2L 79 G 0:0:0:15:0:0 0:0:0:38:0:0 
#	     2L 80 A 12:0:0:0:0:0 38:0:0:0:0:0 
#	     2L 81 A 14:0:0:0:0:0 43:0:0:0:0:0 
#	     2L 82 A 14:0:0:0:0:0 42:0:0:0:0:0
#
#    	     column 1: reference contig
#            column 2: position in the reference contig
#            column 3: refernce character
#            column >3: allele frequencies for all populations in the form A-count:T-count:C-count:G-count:N-count:deletion-count
#
#	     The program accepts the VCF files listed in a CSV files as:
#	     [SAMPLE_NAME],[PATH_TO_VCF_FILE]
#
#	     for example,
#
#	     $ cat vcf_files_MDS42.csv
#	     JC7,../../MDS42/Sample_1-JC7/Sample_1-JC7_lofreq.vcf
#	     JC8,../../MDS42/Sample_2-JC8/Sample_2-JC8_lofreq.vcf
#	     JC9,../../MDS42/Sample_3-JC9/Sample_3-JC9_lofreq.vcf
#	     JC10,../../MDS42/Sample_4-JC10/Sample_4-JC10_lofreq.vcf
#	     JC11,../../MDS42/Sample_5-JC11/Sample_5-JC11_lofreq.vcf
#	     JC12,../../MDS42/Sample_6-JC12/Sample_6-JC12_lofreq.vcf
#		
# Dependencies:
#	     Biopython is installed
#	     PyVCF https://github.com/jamescasbon/PyVCF is installed
#	     Tested on Python 2.6.6
# Version:   0.1 (2018-10-06)
# Authors:   Umer Zeeshan Ijaz (Umer.Ijaz@glasgow.ac.uk)
#                 http://userweb.eng.gla.ac.uk/umer.ijaz
# Created:   2018-10-06
# License:   Copyright (c) 2018 Environmental'Omics Lab, University of Glasgow, UK
#
#            This program is free software: you can redistribute it and/or modify
#            it under the terms of the GNU General Public License as published by
#            the Free Software Foundation, either version 3 of the License, or
#            (at your option) any later version.
#
#            This program is distributed in the hope that it will be useful,
#            but WITHOUT ANY WARRANTY; without even the implied warranty of
#            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#            GNU General Public License for more details.
#
#            You should have received a copy of the GNU General Public License
#            along with this program.  If not, see <http://www.gnu.org/licenses/>.
# **************************************************************/

import sys, getopt
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import vcf
def usage():
	print './lofreq2sync.py -v vcf_files_table.csv'
def main(argv):
	vcf_files_table=''
	flag=0
	try:
		opts, args=getopt.getopt(argv,"hv:",["vcf_files"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			usage()
			sys.exit(2)
		elif opt in ("-v", "--vcf_files"):
			vcf_files_table=arg

        vcf_dict={}
	#Load the sample names and VCF files into dictionary
        if vcf_files_table!='':
                with open(vcf_files_table) as vf:
                        for line in vf:
                                rec=line.rstrip().split(',')
                                vcf_dict[rec[0]]=rec[1]
	ref_pos_dict={}
	sample_pos_allele_string_dict={}
	ref_chromosome=''
	for i in vcf_dict.keys():
		vcf_reader=vcf.Reader(filename=vcf_dict[i])
		for record in vcf_reader:
			#initialise tmp_dict with all zeros
			tmp_dict={}
			tmp_dict['A']=0
			tmp_dict['T']=0
			tmp_dict['G']=0
			tmp_dict['C']=0
			tmp_dict['N']=0
			tmp_dict['*']=0
			CHROM=record.CHROM
			REF=record.REF
			POS=str(record.POS)
			ALT="".join(str(e) for e in record.ALT)
			DP=record.INFO['DP']
			AF=record.INFO['AF']
			#calculate alternate allele frequency count	
			ALTCOUNT=int(round(float(AF)*float(DP)))
			#calculate reference allele frequency count
			REFCOUNT=int(DP)-ALTCOUNT
			
			#assign ALTCOUNT and REFCOUNT to tmp_dict{}
			tmp_dict[ALT]=ALTCOUNT
			tmp_dict[REF]=REFCOUNT
			allele_string=str(tmp_dict['A'])+":"+str(tmp_dict['T'])+":"+str(tmp_dict['G'])+":"+str(tmp_dict['C'])+":"+str(tmp_dict['N'])+":"+str(tmp_dict['*'])	
			#save position strings
			ref_pos_dict[int(POS)]=REF
			#save sample allele_string
			sample_pos_allele_string_dict[i+"#"+POS]=allele_string
			#save reference chromosome name
			ref_chromosome=CHROM

	#now print the actual format
	print "chro\tpos\tref\t"+"\t".join(vcf_dict.keys())
	for i in sorted(ref_pos_dict.keys()):
		sys.stdout.write(ref_chromosome+"\t"+str(i)+"\t"+ref_pos_dict[i])
		for j in vcf_dict.keys():
			sys.stdout.write("\t"+sample_pos_allele_string_dict.get(j+"#"+str(i),"0:0:0:0:0:0"))
		sys.stdout.write("\n")	
if __name__ == "__main__":
	main(sys.argv[1:])

