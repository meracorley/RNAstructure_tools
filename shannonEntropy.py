#script to return the base pairing prob per base given a partition function.
#Adapted from getBPprobs.py from the benchmark project
#Input file format is like: base1	base2	prob_base1_paired_with_base2
# (Input needs to be tab delimited)

import sys
import math
import getopt 


def shannonEntropy(info,size,sqrt=True,log=False): 
    #info is an array containing 3 col for each row: base1 base2 pr(base1-base2)
    #size is the length of the sequence with bpprobs data
    #sqrt and log indicate whther the bpprobs are in sqrt or -log10 format. Can't be both, but can be neither.
    probs = []
    shannonList = []
    #print(info[0:10], len(info), size)
    for i in range(0,size): #creating a nxn matrix to hold base pair probs, where n is seq length
        probs.append([0]*size)
    for i in range(0,len(info)):
        base1 = int(info[i][0])
        base2 = int(info[i][1])
        prob = float(info[i][2])
        if log:
            prob = 10**(-1*float(info[i][2]))
        elif sqrt:
            prob = float(info[i][2])**2.
        probs[base1-1][base2-1] = prob
        probs[base2-1][base1-1] = prob

    for i in range(0, size): #going through every row in the bp probs matrix
        shannon = 0.0 #shannon entropy
        unpaired = 1.0
        for j in range(0, size):
            if(probs[i][j]!=0): #can't take the log of 0
                shannon -= probs[i][j]*math.log(probs[i][j])
                unpaired -= probs[i][j] 
        if(unpaired<=0):
                unpaired = 0.0
        else:
                unpaired = unpaired*math.log(unpaired)
        shannon -= unpaired

        shannonList.append(shannon)
    return shannonList
    

if __name__ == "__main__":
    log = False
    sqrt = False
    inputName = ""
    #############Begin command line options code
    argv = sys.argv[1:] #grabs all the arguments

    def usage():
        print('''python shannonEntropy.py -i input_posteriors.txt [-l || -s] > output.probs
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

    #script to take two strings of dot bracket plots and print out their pearson correlation
    infile = open(inputName, 'r')
    allprobs = infile.read().splitlines() #contains base1 base2 pr_pairing(base1-base2)
    infile.close()
    bpprobs = []
    size = int(allprobs[0].rstrip()) #size of the sequence must be given on the first line of file
    for i in range(2,len(allprobs)): #first two lines of input file should have seq size then a header line.
        thisline = allprobs[i].split()[0:3]
        bpprobs.append(thisline)
    shannon = shannonEntropy(bpprobs,size,sqrt,log)
    for i in shannon:
        print(i)
