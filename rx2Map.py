import numpy as np
import sys

#Converts a .reactivies file to a .map file
#Format: basenum SHAPEval SE Base
#Averages across columns (replicates)
#Also assumes the same number of rows in each column, and that each value is a number.
#USAGE: python rx2SHAPE.py input.txt > output.txt


def getSeqFromFasta(fasta,tnx):
        #Script to extract a given transcript (tnx) from a large fasta file.
        #TO RUN: python getSeqFromFasta.py file.fasta name_of_seq_I_want > output
        #The name searched for shouldn't have spaces in it.
        file = open(fasta, 'r')
        lines = file.read().splitlines()
        file.close()
        i=0 #line counter
        printSeq = False
        seq=""
        found = False #Dont need to keep searching if I already found it.
        while i < len(lines) and not found:
                if lines[i].split()[0][0]==">": #only looking at sequence header lines
                        name=lines[i].split()[0]
                        if name==">"+tnx: #found the sequence, now collect every line
                                printSeq = True
                        elif printSeq: #If I was printing out the sequence I know I already found it. Dont need to search anymore.
                                found = True
                                printSeq = False
                        else:
                                printSeq = False
                elif printSeq:
                        seq+=lines[i].upper()
                
                i+=1            
        return seq

def getAverage(reactivities,name): #return the average SHAPE reactivity of a region
    #print("getAverage")
    seq = getSeqFromFasta("/home/mcorley/annotations/refGene_hg38.v2.fa",name)
    counter = 1
    for i in range(0,len(reactivities)):
        values = reactivities[i,reactivities[i]!=-999]
        if len(values)!=0: #all the values were -999
            average = np.mean(values)
            var = np.var(values)
        else:
            average = -999.
            var = 0
        print(str(counter)+"\t"+str(average)+"\t"+str(var)+"\t"+seq[i])
        counter+=1
if __name__ == "__main__":
    file = open(sys.argv[1], 'r')
    lines = file.read().splitlines()
    file.close()
    if len(lines)==0:
        sys.exit()
    num_columns = len(lines[0].split())
    reactivities = np.zeros((len(lines),num_columns), dtype='float')
    for i in range(0, len(lines)):
        line = lines[i].split()
        for j in range(0,num_columns):
            reactivities[i][j] = float(line[j])
    name = sys.argv[1].split('/')[-1].split('.')[0]
    getAverage(reactivities,name)
