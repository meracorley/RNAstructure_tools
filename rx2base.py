#Converts .reactivity files in the current directory or directory in filepath to a fasta file,
#where each reactivity is converted to a base {A,C,T,G} based on thresholds:
#<0.295953316563; 0.295953316563-0.485405659433; 0.485405659433-0.770302106094, >0.770302106094
#(see plotReactivitiesTnx.ipynb)
#TO RUN: python rx2base [path to rx files (optional)] > output.fa

import numpy as np
import sys
import os

def smooth(profile1d,win):
    if len(profile1d)<win:
        return np.array([])
    smoothed = np.zeros(len(profile1d)-win,dtype='float')
    for i in range(0,len(profile1d)-win):
        thiswin = profile1d[i:i+win][profile1d[i:i+win]!=-999]
        if len(thiswin)==0:
            smoothed[i] = -999.
        else:
            smoothed[i] = np.mean(thiswin)
    return smoothed

def rx2base(rxFile,thresholds):
    file = open(rxFile,'r')
    fileLines = file.read().splitlines()
    file.close()
    start = -1 #I want to know where the first non-999 value occurs and the last
    last = len(fileLines)

    values = np.zeros(len(fileLines),dtype='float')
    for i in range(0,len(fileLines)):
        line = np.array(fileLines[i].split(),dtype='float')
	non999 = line[line!=-999]
        if len(non999)!=0:
            values[i] = np.mean(non999)
            if start==-1:
                start = i
            last = i
        else:
            values[i] = '-999'
    #print(values)
    trimmed = values[start:last+1] #'trim' off ends with no data (keep the middle with good data)
    #print(trimmed)
    smoothed = smooth(trimmed,4) #smooth profile with window size of 4
    numNan = len(smoothed[smoothed==-999])
    if numNan>=5 or len(smoothed)<20:
        return ""
    baseDict = ['A','C','T','G']
    bases = "" #np.zeroes(len(smoothed),dtype='str')    
    for i in range(0,len(smoothed)):
        if smoothed[i]==-999: #if only a few -999 values, replace them with a random pick
            bases+=baseDict[np.random.randint(0,4)]
        else:
            thisbase = 'G'
            notAdded = True
            counter = 0 #ex thresholds = [0.29595,0.4854,0.7703]
            for j in range(0,len(thresholds)):
                if smoothed[i] < thresholds[j]:
                    thisbase = baseDict[j]
                    break
            bases+=thisbase

    return bases

if __name__ == '__main__':
    thresholds = [0.17579,0.42142,0.6571] #these are for smoothed vitro, for non-smoothed vitro: [0.29595,0.4854,0.7703]
    filepath = ""
    if len(sys.argv)>1:
         filepath = sys.argv[1] #can provide the filepath to the desired directory, or current directory is used
    files = [filepath+f for f in os.listdir(filepath+'.') if f.split('.')[-1]=="reactivities" ]
    for file in files:
        sequence = rx2base(file,thresholds)
        name = file.split('/')[-1].split('.')[0]
        if len(sequence)!=0: #rx2base returns empty string for inputs with too many -999 values
            print(">"+name)
            print(sequence)
