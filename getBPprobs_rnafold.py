def main():
	#script to return the base pairing prob per base given a base pairing probability matrix
	#(dp.ps) file from RNAfold -p

	import sys
	import math
	import getopt 
	import os

	sqrt = True
	inputName = ""
	size = 300
	#############Begin command line options code
	argv = sys.argv[1:] #grabs all the arguments

	def usage():
		print('''python getBPprobs.py -i name_dp.ps > name.probs
        Options:
        -h, --help

	REQUIRED
	-i, --input				The _dp.ps file of posterior probabilities from RNAfold.
        ''')

	if len(argv)==0:
		usage()
		sys.exit()	

	try:
		opts, args = getopt.getopt(argv, "hi:", ["help","input="])
	except getopt.GetoptError: #option help input sequence_length
		usage()
		sys.exit(2)

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-i", "--input"):
			inputName = arg

	

	source = "".join(args)
	if len(source)>0:
		print("WARNING: Unused options", source)
	################# END command line options code

	if not os.path.isfile(inputName):
		print(inputName, "does not exist")
		sys.exit()
	infile = open(inputName, 'r')
	allprobs = infile.read().splitlines()
	if len(allprobs)==0:
		print(inputName, "is empty file.")
		sys.exit()
	probs = []
	probStart = -1
	probStop = -1
	for i in range(0,len(allprobs)):
		thisline = allprobs[i].split(' ')
		if len(thisline)>0:
			if thisline[0] == "/sequence":
				size = 0
				i+=1
				while allprobs[i]!=") } def":
					size += len(allprobs[i])-1
					i+=1
		#print("Size is",size)
		if len(thisline)>3 and thisline[3] == "ubox":
			if probStart != -1:
				probStop = i+1
			else:
				probStart = i
	for i in range(0,size): #creating a nxn matrix to hold base pair probs, where n is seq length
		probs.append([0]*size)

	for i in range(probStart,probStop): #starting at line2 because assuming .post file has two header lines
		info = allprobs[i].split()
		base1 = int(info[0])
		base2 = int(info[1])
		prob = float(info[2])**2.		
		probs[base1-1][base2-1] = prob
		probs[base2-1][base1-1] = prob

	#print probs
	#for i in range(0, len(probs[0])):
	#	print(numpy.sum(probs[i]))

	for i in range(0, len(probs)): #going through every row in the bp probs matrix
		downstream = 0.0 #prob of being upstream paired
		for j in range(0, i):
			downstream += probs[i][j]

		upstream = 0.0 #prob of being downstream paired
		for j in range(i+1, size):
			upstream += probs[i][j]

		if upstream > 1.0: #in some cases, rounding errors in UNAfold made bp prob sum > 1.0
			upstream = 1.0
			downstream = 0.0
		elif downstream > 1.0:
			upstream = 0.0
			downstream = 1.0
    
		none = 1.0 - upstream - downstream #prob of not being paired

		#print(str(upstream) + "\t" + str(downstream) + "\t" + str(none)) #comment out if want total pairing prob
		print(upstream+downstream) #comment out if want upstream and downstream probs separately
	infile.close()
main()
