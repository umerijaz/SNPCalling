#!/usr/bin/python
# ***************************************************************
# Name:      gbCompareGenes.py
# Purpose:   Looks at genbank files and returns a presence/absence table for genes
#	     that are present in the list of genbank file provided in a CSV format
#
#	     $ ./gbCompareGenes.py -g gb_files.csv
#		GENE	MDS42	MG1655
#		carB 	1 	1 
#		folA 	1 	1 
#		ksgA 	1 	0 
#		thiP 	1 	1 
#		ddl 	1 	0 
#		aceF 	1 	1 
#		lpxD 	1 	1 
#		metQ 	1 	1 
#		rrf 	1 	0 
#		fadE 	1 	1 
#		frsA 	1 	1 
#		cynX 	1 	1 
#		lacY 	1 	1 
#		lacZ 	1 	1 
#		lacI 	1 	1 
#		thrL |EcoGene:EG11277;GeneID:944742 	0 	1 
#		thrA |EcoGene:EG10998;GeneID:945803 	0 	1 
#		thrB |EcoGene:EG10999;GeneID:947498 	0 	1 
#		thrC |EcoGene:EG11000;GeneID:945198 	0 	1 
#		yaaX |EcoGene:EG14384;GeneID:944747 	0 	1 
#		yaaA |EcoGene:EG10011;GeneID:944749 	0 	1 
#		yaaJ |EcoGene:EG11555;GeneID:944745 	0 	1 
#		talB |EcoGene:EG11556;GeneID:944748 	0 	1 
#		mog |EcoGene:EG11511;GeneID:944760 	0 	1 
#		satP |EcoGene:EG11512;GeneID:944792 	0 	1 
#		yaaW |EcoGene:EG14340;GeneID:944771 	0 	1 
#
#            The program accepts the genbank files listed in a CSV files as:
#            [GENOME_NAME],[PATH_TO_GENBANK_FILE]
#
#            for example,
#
#	     $ cat gb_files.csv
#	     MDS42,../Reference/MDS42reference.gb
#	     MG1655,../Reference/MG1655reference.gb
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
	print './gbCompareGenes.py -g genbank_files_table.csv'

def main(argv):
	genbank_files_table=''
	try:
		opts, args=getopt.getopt(argv,"hg:",["genbank_files"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			usage()
			sys.exit(2)
		elif opt in ("-g", "--genbank_files"):
			genbank_files_table=arg

        gb_dict={}
	#Load the genbank files and their path into the dictionary
        if genbank_files_table!='':
                with open(genbank_files_table) as gf:
                        for line in gf:
                                rec=line.rstrip().split(',')
                                gb_dict[rec[0]]=rec[1]


	unique_genes_id=[]
	unique_genes_details={}
	gb_genes={}


	for i in gb_dict.keys():
		r=SeqIO.read(gb_dict[i],"genbank")
		for feature in r.features:
			if feature.type in ("gene"):
				gene_name="".join(str(e) for e in feature.qualifiers.get("gene",""))
				if len(gene_name)>0:
                                	gene_details=feature.qualifiers.get("db_xref",[])
					if gene_name not in unique_genes_id:
						unique_genes_id.append(gene_name)
						unique_genes_details[gene_name]=gene_details
					gb_genes[i+":"+gene_name]=1

	print "GENE\t"+"\t".join(gb_dict.keys())
	for i in unique_genes_id:
		print i,
		details=unique_genes_details.get(i)
		if details:
			print "|"+";".join(str(e) for e in unique_genes_details.get(i)),
		for j in gb_dict.keys():
			print "\t"+str(gb_genes.get(j+":"+i,0)),
		print "\n",

if __name__ == "__main__":
	main(sys.argv[1:])

