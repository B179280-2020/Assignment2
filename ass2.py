#!/usr/bin/python3
import os, sys, shutil,re
import subprocess
import matplotlib.pyplot as plt


taxon_gp = ""
#ask the user about what protein family and taxonomic group they like
def getInput():
	global taxon_gp
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
seq_Num = findSeq()
print("There are " + str(seq_Num) + " sequences in the dataset you have chosen")
#Test if the the number of the sequences are appropriate and give the user options to decide if they want to continue

def sequence_check(seq_Num):
	while seq_Num <= 1 or seq_Num >= 1000:
		print("Sorry,there is no result or over 1000 sequences. Please check and input again!")
		if seq_Num >= 1 and seq_Num <= 1000:
                        break
		getInput()
		seq_Num = findSeq()
		return seq_Num
sequence_check(seq_Num)

#Tell the user the sequence number and ask them if they are willing to continue
print("There are " + str(seq_number) + " sequences in the dataset you have chosen")
choice1 = input("Do you want to continue?,Y/N\n")
while choice1 == "N":
	if choice1 == "Y":
		break
	print("Thanks, please then input what you want again")
	getInput()
	seq_Num = findSeq()
	print("There are " + str(seq_Num) + " sequences in the dataset you have chosen")
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
spe_N = findSpec()
#Tell the user the number of species in this dataset and ask them if you want to continue
print("There are " + str(spe_N) + " species in this dataset")
choice2 = input("Do you want to continue with the current dataset? Y/N\n")
while choice2 == "N":
        if choice2 == "Y":
                break
        print("Thanks, please then input what you want again")
	getInput()
	seq_Num = findSeq()
	spe_N = findSpec()
        print("There are " + str(spe_N) + " species in this dataset")
        choice2 = input("Do you want to continue?,Y/N\n")


#move on to the main data processing procedure
#using clustalo to align the data
clus = "clustalo -i seq.fa > align.fa"
print(clus)
subprocess.call(clus,shell=True)
print("Clustalo has been successfully done!")
print("The alignment data has been store in the align.fa file")


#using clustalo to align the data
os.system("clustalo -i seq.fa > align.fa")
print("The alignment data has been store in the align.fa file")

#find a representative sequence for blast analysis
#find the sequence with the lease number of "-" in the align.fa file
def _findId():
	f = open('align.fa','r')

	item_list = f.read().split('>')[1:]
	_num = []

	for item in item_list:
		_num.append(item.count('-'))
	print('- number is ' + str(_num[_num.index(min(_num))]) + '\n')
	print('The representative sequence ID is:\n')
	rep_seq = item_list[_num.index(min(_num))]
	ind = rep_seq.find(" ")
        seq_id = rep_seq[0:ind]
        return seq_id

seq_Id=_findId()

#download the test sequence using the seq_id
os.system("esearch -db protein -query " + seq_Id + " | efetch -db protein -format fasta > test_seq.fa")

#using this representative sequence for blast
#make blast data base first
mdb = "makeblastdb -in seq.fa -dbtype prot -out " + taxon_gp
print(mdb)
subprocess.call(mdb,shell=True)


#doing blast using test sequence and data base
blt = "blastp -db " + taxon_gp + " -query test_seq.fa -outfmt 7 > blastoutput.out"
print(blt)
subprocess.call(blt,shell=True)

	
#find the first 250 sequence id of the blastoutput file
fd = "grep -v \"#\" blastoutput.out | cut -f2 | head -250 > blast250.txt"
print(fd)
subprocess.call(fd,shell=True)


#using pullseq to download these 250 sequences
pu = "/localdisk/data/BPSM/Assignment2/pullseq -i seq.fa -n blast250.txt > seq_pull_250.fasta"
print(pu)
subprocess.call(pu,shell=True)

#plotcon to plot the level of conservation
plt = "plotcon -sequence seq_pull_250.fasta -winsize 6 -graph svg"
print(plt)
subprocess.call(plt,shell = True)
print("The plot of the conservation level has been saved as a file.")

#Ask the user if they want to show the plot on the screen
choice3 = input("Do you want to display the plot?,Y/N\n")
if choice3 == 'Y':
    subprocess.call("display plotcon.svg", shell=True)
else:
    print("It's OK, you can check this plot later!")


#Extra EMBOSS analysis to obtain more information
#Using pepstats to obtain more information about the statistics of protein properties
pep = "pepstats -sequence seq_pull_250.fasta -outfile seq250.pepstats -aadata -mwdata -pkdata"
print(pep)
subprocess.call(pep,shell=True)
print("The statistics of protein properties of the 250 sequences from BLAST output have been stored in the seq250.pepstats file.")
print("These statistics include molecular weight, charge, Isoelectric point etc...")



#Using Inforalign analysis 
inf = "infoalign -sequence seq_pull_250.fasta -outfile seq250.infoalign"
print(inf)
subprocess.call(inf,shell = True)
print("The infoalign analysis results have been saved in seq250.infoalign file")


#PROSITE analysis to find motifs
blt250 = open("blast250.txt").read().split("\n")
print(blt250[:-1])
motif_array = []
found_motif = open('fnd_motifs.txt',"w")
found_motif.write("Accession ID\tMotif\n")
for i in blt250[:-1]:
        ese = "esearch -db protein -query " + i + " | efetch -db protein -format fasta > int_seq.fa"
        #print(ese)
        subprocess.call(ese,shell=True)
        pam = "patmatmotifs -sequence int_seq.fa -outfile int_seqmotif -full"
        #print(pam)
        subprocess.call(pam,shell=True)
#find if there are any motifs
        pro_file = open("int_seqmotif").readlines()
        for line in pro_file:           #through each lines
                if re.search('#',line):         #skipping comments
                        next
                elif re.search('Motif', line):       #searching for line where Motif word used
                        line.rstrip()                #cutting new line from the end of line, in order to not to get in my match
                        index = line.find("=")       #motif name is after the sign '=', that's why I want to get index of this sign
                        motif = line[index+2:]       #getting name of the Motif after '= ' and space, so I add 2
                        motif_array.append(motif)    #saving found motif names into array
                        found_motif.write('{0}\t{1}'.format(i,motif))

print("The output of all sequences' scanning for motifs from the PROSITE database are saved in the fnd_motifs.txt file")

