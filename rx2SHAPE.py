import numpy as np
import sys

#Converts a .reactivies file to a .shape file
#Averages across columns (replicates)
#Also assumes the same number of rows in each column, and that each value is a number.
#USAGE: python rx2SHAPE.py input.txt > output.txt

def getAverage(reactivities): #return the average SHAPE reactivity of a region
    #print("getAverage")
    counter = 1
    for i in range(0,len(reactivities)):
        values = reactivities[i,reactivities[i]!=-999]
        if len(values)!=0: #all the values were -999
            average = np.mean(values)
        else:
            average = -999.
        print(str(counter)+"\t"+str(average))
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

    getAverage(reactivities)
