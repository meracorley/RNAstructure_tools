def main():
	#script to return the base pairing prob per base given a base pairing probability matrix.
	#Adapted from getBPprobs.py from the benchmark project
	#Input file format is like: base1	base2	prob_base1_paired_with_base2
	# (Input needs to be tab delimited)

	import sys
	import math
	import getopt 

	log = False
	sqrt = False
	inputName = ""
	
	#############Begin command line options code
	argv = sys.argv[1:] #grabs all the arguments

	def usage():
		print('''python getBPprobs.py -i input_posteriors.txt [-l || -s] > output.probs
        Options:
        -h, --help

	REQUIRED
	-i, --input				The input text file of posterior probabilities. Format is:
						<RNA_size>
						<header>
						base_i	base_j	prob(i paired with j)
						...
        OPTIONAL
        -l, --log                              Probabilities are -log10(prob). RNAstructure uses this format.
	-s, --sqrt				Probabilities are squareRoot(prob). RNAfold uses this format.
        ''')

	if len(argv)==0:
		usage()
		sys.exit()	

	try:
		opts, args = getopt.getopt(argv, "hi:ls", ["help","input=", "log","sqrt"])
	except getopt.GetoptError: #option help input log and sqrt
		usage()
		sys.exit(2)

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-i", "--input"):
			inputName = arg

		elif opt in ("-l", "--log"):
			log=True
		
		elif opt in ("-s", "--sqrt"):
			sqrt=True
	
	if log and sqrt:
		print("Probabilities cannot be in both -log10 and sqrt format.")
		usage()
		sys.exit()

	source = "".join(args)
	if len(source)>0:
		print("WARNING: Unused options", source)
	################# END command line options code

	infile = open(inputName, 'r')
	allprobs = infile.read().splitlines()
	probs = []
	size = int(allprobs[0].rstrip()) #size of the sequence must be given
	for i in range(0,size): #creating a nxn matrix to hold base pair probs, where n is seq length
		probs.append([0]*size)
	for i in range(2,len(allprobs)): #starting at line2 because assuming .post file has two header lines
		info = allprobs[i].split()
		base1 = int(info[0])
		base2 = int(info[1])
		prob = float(info[2])
		if log:
			prob = 10**(-1*float(info[2]))
		elif sqrt:
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

		print(str(upstream) + "\t" + str(downstream) + "\t" + str(none)) #comment out if want total pairing prob
		#print(upstream+downstream) #comment out if want upstream and downstream probs separately
	infile.close()
main()
