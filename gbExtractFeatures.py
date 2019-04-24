#!/usr/bin/python
# ***************************************************************
# Name:      gbExtractFeatures.py
# Purpose:   After SNP calling is done using lofreq, the program can take multiple VCF files
#	     provided in a CSV format and 
#		a) generate a tab delimited list to show which genes are affected by SNPs
#		b) annotate a genbank file with SNPs and produce secondary genbank file(s) (either single or separate) with annotated SNPs
#
#	     At the bare minimal, you can first explore the SNPs using -f 0 switch:
#
#	     $ ./gbExtractFeatures.py -g ../Reference/MDS42reference.gb -v vcf_files_MDS42.csv -f 0 
#		SAMPLE=JC7,  CHROM=NC_020518.1,  POS=61609,  REF=A,  ALT=T,  INFO={'SB': 12, 'DP4': [2337, 2323, 34, 53], 'DP': 4765, 'AF': 0.018258} 
#		SAMPLE=JC7,  CHROM=NC_020518.1,  POS=61621,  REF=T,  ALT=C,  INFO={'SB': 4, 'DP4': [2369, 2375, 36, 45], 'DP': 4834, 'AF': 0.016756} 
#		SAMPLE=JC7,  CHROM=NC_020518.1,  POS=120642,  REF=A,  ALT=G,  INFO={'SB': 56, 'DP4': [2018, 2272, 84, 38], 'DP': 4414, 'AF': 0.027639} ,  
#				TYPE=gene,  PRODUCT=None,  LOCUS_TAG=ECMDS42_RS00555,  PROTEIN_ID=None,  TRANSLATION=None ,  TYPE=CDS,  
#				PRODUCT=['dihydrolipoyllysine-residue acetyltransferase component of pyruvate dehydrogenase complex'],  
#				LOCUS_TAG=ECMDS42_RS00555,  PROTEIN_ID=['WP_000963518.1'],  TRANSLATION=['MAIEIKVPDIGADEVEITEILVKVGDKV
#				EAEQSLITVEGDKASMEVPSPQAGIVKEIKVSVGDKTQTGALIMIFDSADGAADAAPAQAEEKKEAAPAAAPAAAAAKDVNVPDIGSDEVEVTEILVKVGDK
#				VEAEQSLITVEGDKASMEVPAPFAGTVKEIKVNVGDKVSTGSLIMVFEVAGEAGAAAPAAKQEAAPAAAPAPAAGVKEVNVPDIGGDEVEVTEVMVKVGDKV
#				AAEQSLITVEGDKASMEVPAPFAGVVKELKVNVGDKVKTGSLIMIFEVEGAAPAAAPAKQEAAAPAPAAKAEAPAAAPAAKAEGKSEFAENDAYVHATPLIR
#				RLAREFGVNLAKVKGTGRKGRILREDVQAYVKEAIKRAEAAPAATGGGIPGMLPWPKVDFSKFGEIEEVELGRIQKISGANLSRNWVMIPHVTHFDKTDITE
#				LEAFRKQQNEEAAKRKLDVKITPVVFIMKAVAAALEQMPRFNSSLSEDGQRLTLKKYINIGVAVDTPNGLVVPVFKDVNKKGIIELSRELMTISKKARDGKL
#				TAGEMQGGCFTISSIGGLGTTHFAPIVNAPEVAILGVSKSAMEPVWNGKEFVPRLMLPISLSFDHRVIDGADGARFITIINNTLSDIRRLVM'] 
#		SAMPLE=JC7,  CHROM=NC_020518.1,  POS=120645,  REF=C,  ALT=T,  INFO={'SB': 24, 'DP4': [2027, 2248, 64, 39], 'DP': 4392, 'AF': 0.023452000000000001} ,  
#				TYPE=gene,  PRODUCT=None,  LOCUS_TAG=ECMDS42_RS00555,  PROTEIN_ID=None,  TRANSLATION=None ,  TYPE=CDS,  
#				PRODUCT=['dihydrolipoyllysine-residue acetyltransferase component of pyruvate dehydrogenase complex'],  
#				LOCUS_TAG=ECMDS42_RS00555,  PROTEIN_ID=['WP_000963518.1'],  TRANSLATION=['MAIEIKVPDIGADEVEITEILVKVGDKV
#				EAEQSLITVEGDKASMEVPSPQAGIVKEIKVSVGDKTQTGALIMIFDSADGAADAAPAQAEEKKEAAPAAAPAAAAAKDVNVPDIGSDEVEVTEILVKVGDK
#				VEAEQSLITVEGDKASMEVPAPFAGTVKEIKVNVGDKVSTGSLIMVFEVAGEAGAAAPAAKQEAAPAAAPAPAAGVKEVNVPDIGGDEVEVTEVMVKVGDKV
#				AAEQSLITVEGDKASMEVPAPFAGVVKELKVNVGDKVKTGSLIMIFEVEGAAPAAAPAKQEAAAPAPAAKAEAPAAAPAAKAEGKSEFAENDAYVHATPLIR
#				RLAREFGVNLAKVKGTGRKGRILREDVQAYVKEAIKRAEAAPAATGGGIPGMLPWPKVDFSKFGEIEEVELGRIQKISGANLSRNWVMIPHVTHFDKTDITE
#				LEAFRKQQNEEAAKRKLDVKITPVVFIMKAVAAALEQMPRFNSSLSEDGQRLTLKKYINIGVAVDTPNGLVVPVFKDVNKKGIIELSRELMTISKKARDGKL
#				TAGEMQGGCFTISSIGGLGTTHFAPIVNAPEVAILGVSKSAMEPVWNGKEFVPRLMLPISLSFDHRVIDGADGARFITIINNTLSDIRRLVM'] 
#
#            The program accepts the genbank files listed in a CSV files as:
#            [SAMPLE_NAME],[PATH_TO_VCF_FILE]
#
#            for example,
#
#	     $ cat vcf_files_MDS42.csv
#	       JC7,../MDS42/Sample_1-JC7/Sample_1-JC7_lofreq.vcf
#	       JC8,../MDS42/Sample_2-JC8/Sample_2-JC8_lofreq.vcf
#	       JC9,../MDS42/Sample_3-JC9/Sample_3-JC9_lofreq.vcf
#	       JC10,../MDS42/Sample_4-JC10/Sample_4-JC10_lofreq.vcf
#	       JC11,../MDS42/Sample_5-JC11/Sample_5-JC11_lofreq.vcf
#	       JC12,../MDS42/Sample_6-JC12/Sample_6-JC12_lofreq.vcf
#               
# Dependencies:
#            Biopython is installed
#            PyVCF https://github.com/jamescasbon/PyVCF is installed
#            Tested on Python 2.6.6
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
	print './gbExtractFeatures.py -g genbank_file.gb -v vcf_files_table.csv -f FLAG'
	print '-f 0(default)	print features'
	print '-f 1		annotate and save genebank files separately'
        print '-f 2             annotate and save compound genbank file with unique features'
        print '-f 3             annotate and save compound genbank file with all features'
def annotate_genbank_file_with_features(genbank_file,vcf_dict):
        for i in vcf_dict.keys():
                vcf_reader=vcf.Reader(filename=vcf_dict[i])
                r=SeqIO.read(genbank_file,"genbank")
                for record in vcf_reader:
                        my_feature_location=FeatureLocation(record.POS,record.POS,strand=1)
                        my_feature_type="SNP"
                        my_feature_qualifiers={"ALT":''.join(str(e) for e in record.ALT),"AF":str(record.INFO['AF'])}
                        my_feature=SeqFeature(my_feature_location,type=my_feature_type,qualifiers=my_feature_qualifiers)
                        r.features.append(my_feature)
                SeqIO.write(r,i+".gb","genbank")	
def annotate_compound_genbank_file_unique_features(genbank_file,vcf_dict):
	r=SeqIO.read(genbank_file,"genbank")
	unique_positions=[]
	for i in vcf_dict.keys():
		vcf_reader=vcf.Reader(filename=vcf_dict[i])
                for record in vcf_reader:
                        my_feature_location=FeatureLocation(record.POS,record.POS,strand=1)
                        my_feature_type="SNP"
                        my_feature_qualifiers={"ALT":''.join(str(e) for e in record.ALT),"AF":str(record.INFO['AF'])}
                        my_feature=SeqFeature(my_feature_location,type=my_feature_type,qualifiers=my_feature_qualifiers)
			if record.POS not in unique_positions:
                        	r.features.append(my_feature)
				unique_positions.append(record.POS)
	SeqIO.write(r,"compound_unique.gb","genbank")
def annotate_compound_genbank_file_all_features(genbank_file,vcf_dict):
        r=SeqIO.read(genbank_file,"genbank")
        for i in vcf_dict.keys():
                vcf_reader=vcf.Reader(filename=vcf_dict[i])
                for record in vcf_reader:
                        my_feature_location=FeatureLocation(record.POS,record.POS,strand=1)
                        my_feature_type="SNP"+"_"+i
                        my_feature_qualifiers={"ALT":''.join(str(e) for e in record.ALT),"AF":str(record.INFO['AF'])}
                        my_feature=SeqFeature(my_feature_location,type=my_feature_type,qualifiers=my_feature_qualifiers)
                        r.features.append(my_feature)
        SeqIO.write(r,"compound_all.gb","genbank")
def print_features(genbank_file,vcf_dict):
        for r in SeqIO.parse(genbank_file,"genbank"):
                for i in vcf_dict.keys():
                        vcf_reader=vcf.Reader(filename=vcf_dict[i])
                        for record in vcf_reader:
                                print "SAMPLE="+i+", ",
                                print "CHROM="+record.CHROM+", ",
                                print "POS="+str(record.POS)+", ",
                                print "REF="+''.join(record.REF)+", ",
                                print "ALT="+''.join(str(e) for e in record.ALT)+", ",
                                print "INFO="+str(record.INFO),
                                for feature in r.features:
                                        if record.POS in feature:
                                                if feature.type in ("CDS","gene"):
                                                        print ", ",
                                                        print "TYPE="+str(feature.type)+", ",
                                                        print "PRODUCT="+str(feature.qualifiers.get("product"))+", ",
                                                        print "LOCUS_TAG="+str("".join(feature.qualifiers.get("locus_tag")))+", ",
                                                        print "PROTEIN_ID="+str(feature.qualifiers.get("protein_id"))+", ",
                                                        print "TRANSLATION="+str(feature.qualifiers.get("translation")),

                                print "\n",


def main(argv):
	genbank_file=''
	vcf_files_table=''
	flag=0
	try:
		opts, args=getopt.getopt(argv,"hg:v:f:",["genbank","vcf_files","flag"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			usage()
			sys.exit(2)
		elif opt in ("-g", "--genbank"):
			genbank_file=arg
		elif opt in ("-v", "--vcf_files"):
			vcf_files_table=arg
		elif opt in ("-f", "--flag"):
			flag=int(arg)

        vcf_dict={}
	#Load the sample names and VCF files into dictionary
        if vcf_files_table!='':
                with open(vcf_files_table) as vf:
                        for line in vf:
                                rec=line.rstrip().split(',')
                                vcf_dict[rec[0]]=rec[1]
	if flag==0:
		print_features(genbank_file,vcf_dict)					
	elif flag==1:
		annotate_genbank_file_with_features(genbank_file,vcf_dict)
	elif flag==2:
		annotate_compound_genbank_file_unique_features(genbank_file,vcf_dict)
	elif flag==3:
		annotate_compound_genbank_file_all_features(genbank_file,vcf_dict)
if __name__ == "__main__":
	main(sys.argv[1:])

