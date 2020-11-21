#!/usr/bin/python3
import os, sys, shutil
import subprocess
import matplotlib.pyplot as plt

#ask the user about what protein family and taxonomic group they like
def getInput():
	protein_fam = input("Enter the protein family name\n")
	taxon_gp = input("Enter the taxonomic group name\n")
	print("Thanks, you have chosen " + protein_fam + " in " + taxon_gp + "\n")
	es1 = 'esearch -db protein -query \" '+ taxon_gp +'[orgn] AND '+ protein_fam + '[title]\" | efetch -db protein -format fasta > seq.fa '
	print("This is what will be run: " + es1)
#Run the esearch and efetch command to get the dataset from NCBI
	return subprocess.call(es1,shell=True)
getInput()
print("The protein sequences of this taxonomic group have been stored in the seq.fa")

#Find how many sequences are there in the dataset user have chosen
def findSeq():
	seqnumber= subprocess.getoutput("grep -c \">\" seq.fa")
	seq_number = int(seqnumber)
	return seq_number
findSeq()
print("There are " + str(seq_number) + " sequences in the dataset you have chosen")
#Test if the the number of the sequences are appropriate and give the user options to decide if they want to continue

def sequence_check(seq_number):
	while seq_number <= 1 or seq_number >= 1000:
		print("Sorry,there is no result or over 1000 sequences. Please input again!")
		if seq_number >= 1 and seq_number <= 1000:
                        break
		getInput()
		findSeq()
		return seq_number
sequence_check(seq_number)

#Tell the user the sequence number and ask them if they are willing to continue
print("There are " + str(seq_number) + " sequences in the dataset you have chosen")
choice1 = input("Do you want to continue?,Y/N\n")
while choice1 == "N":
	if choice1 == "Y":
		break
	print("Thanks, please then input what you want again")
	getInput()
	findSeq()
	print("There are " + str(seq_number) + " sequences in the dataset you have chosen")
	choice1 = input("Do you want to continue?,Y/N\n")
	
#find the number of species in this dataset
def findSpec():
	get_header = "grep \">\" seq.fa > seq_header.fa"
	subprocess.call(get_header,shell=True)
	seq = open("seq_header.fa")
	spe_num = []
	for i in seq:
		l1 = i.find('[')
		l2 = i.find(']')
		spe_num.append(i[l1:l2])
	spe_n = len(set(spe_num))
	return spe_n
findSpec()
#Tell the user the number of species in this dataset and ask them if you want to continue
print("There are " + str(spe_n) + " species in this dataset")
choice2 = input("Do you want to continue with the current dataset? Y/N\n")
while choice2 == "N":
        if choice2 == "Y":
                break
        print("Thanks, please then input what you want again")
        getInput()
        findSeq()
	findSpec()
        print("There are " + str(spe_n) + " species in this dataset")
        choice2 = input("Do you want to continue?,Y/N\n")

#move on to the main data processing procedure

#using clustalo to align the data
os.system("clustalo -i seq.fa > align.fa")
print("The alignment data has been store in the align.fa file")

#find a representative sequence for blast analysis
def _find():
	f = open('align.fa','r')

	item_list = f.read().split('>')[1:]
	_num = []

	for item in item_list:
		_num.append(item.count('-'))

	print('- number is ' + str(_num[_num.index(min(_num))]) + '\n')
	print('The representative sequence is:\n')

	return item_list[_num.index(min(_num))]	

_find()
os.system("esearch -db protein -query seq_id | efetch -db protein -format fasta > test_seq.fa")
#using this representative sequence for blast
#make blast data base first
subprocess.call("makeblstdb -in seq.fa -dbtype prot -out " + taxon_gp)

#doing blast using test sequence and data base


	



