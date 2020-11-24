#!/usr/bin/python3
import os, sys, shutil,re
import subprocess

taxon_gp = ""
#ask the user about what protein family and taxonomic group they like
def getInput():
	global taxon_gp           #set it as a global varibale
	protein_fam = input("Enter the protein family name\n")            #ask for an input
	taxon_gp = input("Enter the taxonomic group name\n")
	print("Thanks, you have chosen " + protein_fam + " in " + taxon_gp + "\n")
	es1 = 'esearch -db protein -query \" '+ taxon_gp +'[orgn] AND '+ protein_fam + '[title]\" | efetch -db protein -format fasta > seq.fa '
	print("This is what will be run: " + es1)
#Run the esearch and efetch command to get the dataset from NCBI
	return subprocess.call(es1,shell=True)
getInput()                   #run this function
print("The protein sequences of this taxonomic group have been stored in the seq.fa")

#Find how many sequences are there in the dataset user have chosen
def findSeq():
	seqnumber= subprocess.getoutput("grep -c \">\" seq.fa")              #obtain the number of sequences in the downloaded dataset
	seq_number = int(seqnumber)                                          #change it into an integer
	return seq_number
seq_Num = findSeq()                                                          #run this function


#Test if the the number of the sequences are appropriate and give the user options to decide if they want to continue
def sequence_check(seq_Num):
	while seq_Num <= 1 or seq_Num >= 1000:                              #the number of sequences should not be less than 1 or more than 1000
		print("Sorry,there is something wrong with your input - no result or over 1000 sequence. Please input again!")
		if seq_Num >= 1 and seq_Num <= 1000:                        #break the loop if the sequences are in normal range
                        break
		getInput()                                                  #if not in normal range, run the function above to ask for input again
		seq_Num = findSeq()                                         
		return seq_Num
sequence_check(seq_Num)

#Tell the user the sequence number and ask them if they are willing to continue
print("There are " + str(seq_Num) + " sequences in the dataset you have chosen")      #tell the user the seqs number
choice1 = input("Do you want to continue?,Y/N\n")                          #ask the user if they are happy with the current dataset 
while choice1 == "N":                                                      #loop until the user is happy
	if choice1 == "Y":
		break
	print("Thanks, please then input what you want again")
	getInput()
	seq_Num = findSeq()
	sequence_check(seq_Num)
	print("There are " + str(seq_Num) + " sequences in the dataset you have chosen")
	choice1 = input("Do you want to continue?,Y/N\n")
	
#find the number of species in this dataset
def findSpec():
	get_header = "grep \">\" seq.fa > seq_header.fa"        #obtain a new file which contains all the seqences header
	subprocess.call(get_header,shell=True) 
	seq = open("seq_header.fa")
	spe_num = []                                            
	for i in seq:
		l1 = i.find('[')                                #the species' names are in [ ]
		l2 = i.find(']')
		spe_num.append(i[l1:l2])                        #get a list of all species'names
	spe_n = len(set(spe_num))                               #get rid of duplicates
	return spe_n
spe_N = findSpec()                                              #run this function and return species number


#Tell the user the number of species in this dataset and ask them if you want to continue
print("There are " + str(spe_N) + " species in this dataset")
choice2 = input("Do you want to continue with the current dataset? Y/N\n")          #ask the users if they are happy
while choice2 == "N":                                                            #loop until the user is happy
	if choice2 == "Y":
		break
	print("Thanks, please then input what you want again")
	getInput()
	seq_Num = findSeq()
	spe_N = findSpec()
	sequence_check(seq_Num)
	print("There are " + str(seq_Num) + " sequences in the dataset you have chosen")
	print("There are " + str(spe_N) + " species in this dataset")
	choice2 = input("Do you want to continue?,Y/N\n")


#Ask the user if they want to remove partial seqs in the dataset for further analysis
print("Just a reminder, there may have some partial sequences in this dataset.")
choice3 = input("Do you want to remove the partial sequences for further analysis? Y/N\n")        #ask users if they want to remove partial seqs
if choice3 == "Y":
	fu = "grep -v \"partial\" seq_header.fa > seq_fheader.fa"                           #remove partial seqs and put the rest headers into a new file
	#print(fu)
	subprocess.call(fu,shell=True)

#obtain headers without partial seqs
	accf = open("accf.txt","w")
	for i in open("seq_fheader.fa"):
		h = i.split(">")                                            #remove ">" in order to successfully use pullseq
		accf.write(h[1])
	accf.close()	

#pullseq all seqs from the header file generated above
	fl = "/localdisk/data/BPSM/Assignment2/pullseq -i seq.fa -n accf.txt > seq_full.fa"
	#print(fl)
	subprocess.call(fl,shell=True)
	os.remove("seq.fa")
	os.rename("seq_full.fa","seq.fa")                                  #replace the order seq.fa with the new one
	print("Partial sequences have been removed, new sequences dataset has been stored in the seq.fa file")

else:
	print("OK, Let's do some analysis for the current dataset!")

#move on to the main data processing procedure
#using clustalo to align the data
print("Let's do some analysis!")
print("First, Clustalo will be run!")
clus = "clustalo -i seq.fa > align.fa"
#print(clus)
subprocess.call(clus,shell=True)                      #do clustalo
print("Clustalo has been successfully done!")
print("The alignment data has been store in the align.fa file")



#find a representative sequence for blast analysis
#find the sequence with the lease number of "-" in the align.fa file
print("Let's find the representative sequence from the aligment results!")          #BLAST needs a test sequence, thus we need to find a representative one
def _findId():
	f = open('align.fa','r')

	item_list = f.read().split('>')[1:]      #split the contents and make it into a list
	_num = []

	for item in item_list:
		_num.append(item.count('-'))     #count the number of "-" for each seq
	rep_seq = item_list[_num.index(min(_num))]        #find the rep_seq which has the least number of "-"
	ind = rep_seq.find(" ")
	seq_id = rep_seq[0:ind]                  #obtain the representative seq ID
	print('The representative sequence ID is: ' + seq_id)
	return seq_id

seq_Id=_findId()

#download the test sequence using the seq_id
print("Second,BLAST will be run!")
print("The representative sequence will be downloaded and served as a test sequence for BLAST.")
os.system("esearch -db protein -query " + seq_Id + " | efetch -db protein -format fasta > test_seq.fa")    #download the representative sequence

#using this representative sequence for blast analysis
#make blast data base first

mdb = "makeblastdb -in seq.fa -dbtype prot -out " + taxon_gp              #make BLAST database
#print(mdb)
subprocess.call(mdb,shell=True)


#doing blast using test sequence and data base
blt = "blastp -db " + taxon_gp + " -query test_seq.fa -outfmt 7 > blastoutput.out"
#print(blt)
subprocess.call(blt,shell=True)
print("BLAST analysis has been successfully done!")
	
#find the first 250 sequence id of the blastoutput file
print("The 250 most similar sequences will be obtained from the BLAST results!")
fd = "grep -v \"#\" blastoutput.out | cut -f2 | head -250 > blast250.txt"  #based on BLAST output, find the 250 most similar seqs, put the headers in a file
#print(fd)
subprocess.call(fd,shell=True)


#using pullseq to download these 250 sequences
print("Now, it's time to download the 250 most similar sequences using pullseq!")
pu = "/localdisk/data/BPSM/Assignment2/pullseq -i seq.fa -n blast250.txt > seq_pull_250.fasta"       #pullseq can download seqs if input a seq header file
#print(pu)
subprocess.call(pu,shell=True)               #download the 250 most similar seqs

#plotcon to plot the level of conservation
print("Let's plot the level of conservation based on this dataset!")
plt = "plotcon -sequence seq_pull_250.fasta -winsize 6 -graph svg"      #plotcon can make a picture of the level of protein sequence conservation
#print(plt)
subprocess.call(plt,shell = True)
print("The plot of the conservation level has been saved as a file.")

#Ask the user if they want to show the plot on the screen
choice4 = input("Do you want to display the plot?,Y/N\n")          #ask the users if they want to see the plot right now
if choice4 == 'Y':
    subprocess.call("display plotcon.svg", shell=True)             #display the plot if they want to 
else:
    print("It's OK, you can check this plot later!")



#PROSITE analysis to find motifs                    #PROSITE can only accept one sequence per time
blt250 = open("blast250.txt").read().split("\n")          #process the 250 most similar seqs header, we only need the ID of the seq in this step
#print(blt250[:-1])
motif_list = []                                           #make a list for storing motifs
#open a new txt file and make it writable
found_motif = open('fnd_motifs.txt',"w")                 #creat a new file for storing the motifs information, make it writable
found_motif.write("Accession ID\tMotif\n")               #write in title
for i in blt250[:-1]:                                    #the last element of blt250 is "", thus we don't need it
        ese = "esearch -db protein -query " + i + " | efetch -db protein -format fasta > int_seq.fa"         #download each seq based on its ID
        #print(ese)
        subprocess.call(ese,shell=True)
        pam = "patmatmotifs -sequence int_seq.fa -outfile int_seqmotif -full"                      #scan the protein seq with motifs from PROSITE database
        #print(pam)
        subprocess.call(pam,shell=True)
#find if there are any motifs
        pro_file = open("int_seqmotif").readlines()                              #open the output file after patmat
        for line in pro_file:           #go through each line
                if re.search('#',line):         #skipping the first few comments lines
                        next
                elif re.search('Motif', line):       #searching for "Motif" word
                        line.rstrip()               
                        index = line.find("=")       #locate the motif 
                        motif = line[index+2:]       #getting name of the Motif after '= ' and space
                        motif_list.append(motif)    #saving found motif names into the list
                        found_motif.write('{0}\t{1}'.format(i,motif))
found_motif.close()
print("The output of all sequences' scanning for motifs from the PROSITE database are saved in the fnd_motifs.txt file")
choice5 = input("Do you want to see the motif right now? Y/N\n")
if choice5 == "Y":
	print("\n".join(set(motif_list)))
else:
	print("OK, you can check the detailed information in the fnd_motifs.txt file!")


#Extra EMBOSS analysis to obtain more information
#Using pepstats to obtain more information about the statistics of protein properties
print("Let's do pepstats to obtain a more detailed statistics of the proteins")
pep = "pepstats -sequence seq_pull_250.fasta -outfile seq250.pepstats"             #pepstats can obtain some basic statistics of the protein seqs
#print(pep)
subprocess.call(pep,shell=True)
print("The statistics of protein properties of the 250 sequences from BLAST output have been stored in the seq250.pepstats file.")
print("These statistics include molecular weight, charge, Isoelectric point etc...")


#Using Inforalign 
print("Let's do infoalign to get extra information")            
inf = "infoalign -sequence seq_pull_250.fasta -outfile seq250.infoalign"
#print(inf)
subprocess.call(inf,shell = True)
print("The infoalign analysis results have been saved in seq250.infoalign file")
print("The summary information from pepstats and infoalign has been saved, you can check them later!")

#Ending
print("Thanks for using this programme! This programme has finished! Hope you enjoy it!")
