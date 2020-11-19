#!/usr/bin/python3
import os, sys, shutil
import subprocess
import matplotlib.pyplot as plt

#ask the user about what protein family and taxonomic group they like
protein_fam = input("Enter the protein family name\n")
taxon_gp = input("Enter the taxonomic group name\n")
print("Thanks, you have chosen " + protein_fam + " in " + taxon_gp + "\n")
es1 = 'esearch -db protein -query \" '+ taxon_gp +'[organism] AND '+ protein_fam + '[protein]\" | efetch -db protein -format fasta > seq.fa '
print("This is what will be run: " + es1)
subprocess.call(es1,shell=True)
print("The protein sequences of this taxonomic group have been stored in the seq.fa")
f_seq = "grep -c \">\" seq.fa"
seq_num = subprocess.call(f_seq, shell=True)
print("There are " + str(seq_num) + " sequences in the dataset you have chosen.")
get_header = "grep \">\" seq.fa > seq_header.fa"
subprocess.call(get_header,shell=True)
seq = open("seq_header.fa")
spe_num = []
for i in seq:
	l1 = i.find('[')
	l2 = i.find(']')
	spe_num.append(i[l1:l2])
spe_n = len(set(spe_num))
print("There are " + str(spe_n) + " species in this dataset")

	

	



