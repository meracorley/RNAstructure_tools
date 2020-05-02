
#Adapted from profileCorrcoef.py in ~/alain/combined/
#Gives the correlation coefficient between two profiles in windows, rather than using the whole profile.
#Have to have atleast two profiles to compare (duh). Can compare more than 2, but this is only useful in a region that
#all the isoforms have and where there aren't a lot of -999s.
#The profile can be of anything--shape values, BGrates, etc.
#Usage: python windowCorrcoef.py INPUT winsize[optional] > outfile
#Outfile format:
import numpy
import sys
    
def preprocess(datalist): #datalist should be a 2-d array
	profiles = []
	position = [-1]*len(datalist) #for each base in given file will hold its position in corrcoefs, or "-1"
					#For instance, if the input has -999 in the 2nd position, then position=[0,-1,1,2...]
	counter = -1 #to keep track of base positions after -999s are deleted
    
	#This formats each sample to remain as a column; each sample position is a row
	for i in range(0,len(datalist)):
		counter+=1
		thisline = datalist[i]
		numbers = []
		add=True
		for j in thisline:
			value = -999.
			if j!="NA" and j!="":
				value= float(j)
			if value==-999.0:
				position[i]=-1
				counter-=1
				add=False #dont add this row to profiles if one of the samples has a -999
				break
			numbers.append(value)
		if add:
			profiles.append(numbers)
			position[i]=counter
 
	''' #This formats each sample column to be a row in profiles, instead of remaining as a column
	for j in range(0,len(datalist[0].split())):
		profiles.append([])
        
	for i in range(0,len(datalist)): 
		counter+=1
		for j in range(0,len(datalist[0].split())):
			#print datalist[i].split()[j]
			value = datalist[i].split()[j]
			if value == "NA" or value =="":
				value = -999
			value = float(value)
			if value == -999:
				position[i]=-1
				counter-=1
				for k in range(0,j): #need to delete values in the other profiles if this profile has a -999
					del profiles[k][-1]
				break #don't use any of the isoform's values if one of them is -999
			if value <0: #why not keep negative values??
				value = 0.0
			profiles[j].append(value)
			position[i] = int(counter)
	'''
	#print(profiles[-50:-1])
	return profiles, position

def winCorrcoef(profiles, winsize=50, stepsize=1): #the format of profiles should be a 2D array, one sample per column.
	#print profiles
	#print position
	#corrcoef= numpy.corrcoef(profiles)					#number of rows is number of bases (non -999)
	n = len(profiles[0]) #the number of samples, if each sample is in a column of profiles
	#print("n is ", n)
	x = len(profiles) #the number of positions in each sample
	corrcoefs = [[numpy.nan]*int((n)*(n-1)/2) for i in range(x)] #for storing all window-calculated corrcoefs for every profile comparison
	#corrcoefs = [[numpy.nan]*int((n)*(n-1)/2)) for i in range((x-winsize+step)/step)]
	#The number of comparisons is n choose 2, which = (n)(n-1)/2, where len(profiles[0]) gives me n.
	#print("x is", x)
	#print(x-winsize+1-x%stepsize)
	#print(profiles[-50:-1])
	for i in range(0,x-winsize+1-x%stepsize,stepsize): #fyi: total number of corrcoefs taken = floor_(x-winsize+step)/step
		winprofile = []
		for j in range(i,i+winsize): #For each profile/sample, append this window of values to winprofile
			winprofile.append(profiles[j])
			#print "winprofile is", winprofile
		wincorrcoef = numpy.corrcoef(winprofile, rowvar=False) #rowvar=F indicates to compare columns
		#print "wincorrcoef is", wincorrcoef
		combination = 0
		for row in range(0,n): #iterating through the corrcoefs in the top right above the diagonal
			for col in range(row+1,n):
				corrcoefs[i][combination] = wincorrcoef[row][col] #corrcoef between the first and second profiles
				#print "corrcoefs is", corrcoefs
				combination += 1
	return corrcoefs

def returnCorrcoefs(corrcoefs,position): #return corrcoef as an array (instead of printing)
	n = len(corrcoefs[0]) #the number of pairwise corrcoefs
	fullProfile = numpy.zeros((len(position),n),dtype="float")

	corrcoefs.append([numpy.nan]*n) #shortcut to otputting "NA" for positions with -999s that will have "-1" in position list. -1 will refer to this last posi$
	for i in range(0,len(position)): #position tells me the position of each base's corrcoeff value in corrcoefs
		buildstr = []
		for j in range(n):
			buildstr.append(corrcoefs[position[i]][j])
		fullProfile[i] = buildstr
	return fullProfile

#This part prints out the corrcoef profiles#########
def printCorrcoefs(corrcoefs, position):
	n = len(corrcoefs[0]) #the number of pairwise corrcoefs
	corrcoefs.append([numpy.nan]*n) #shortcut to otputting "NA" for positions with -999s that will have "-1" in position list. -1 will refer to this last position appended to corrcoefs.
	for i in position: #position tells me the position of each base's corrcoeff value in corrcoefs
		buildstr = ""
		for j in range(n):
			buildstr += str(corrcoefs[i][j])+"\t"
		print(buildstr.rstrip('\t')) #i itself is an index of a base's corrcoef value in corrcoefs

	'''
	#Calculating the average corrcoef
	for j in range(0, len(corrcoefs[0])): #the different pairwise profile comparison
		sum = 0.
		count = 0.
		for i in range(0, len(corrcoefs)): #the non NA corrcoef positions
			if corrcoefs[i][j]!="NA":
				sum+=corrcoefs[i][j]
				count+=1
		print(sum/count)
	'''     
if __name__ == "__main__":
	data = open(sys.argv[1], 'r') #Should be in format: profile1    profile2...profilen
	lines = data.read().splitlines()
	datalist = numpy.zeros((len(lines),len(lines[0].split())), dtype='object')
	for i in range(0,len(lines)):
		line = lines[i].split()
		datalist[i] = line
	profiles, position = preprocess(datalist)
	winsize = 50 #default
	stepsize = 1
	if len(sys.argv)>2:
		winsize = int(sys.argv[2])
	if len(profiles) < winsize:
		winsize = len(profiles)
	#print winsize
	stepsize = 1
	corrcoefs = winCorrcoef(profiles, winsize, stepsize)
	fullProfile = returnCorrcoefs(corrcoefs, position)
	printCorrcoefs(corrcoefs, position)
	data.close()
