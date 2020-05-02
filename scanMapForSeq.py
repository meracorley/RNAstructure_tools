import sys, os, getopt
def scanForSeq(inputName,pattern,startPos=0,stopPos=0,numMatches=-1):
    if os.path.isfile(inputName): #input is a MaP file (from probing experiments), format is 4columns: basenum rx_value SE nucleotide
        infile = open(inputName, 'r')
        mapfile = infile.read().splitlines()
        infile.close()
    else:
        print(inputName, "does not exist.")
        return []
    seq = ""
    for line in mapfile:
        thisline = line.split()
        if len(thisline)<4:
            print(inputName, "does not have 4 columns at every line.")
            return []
        seq+=thisline[3]  
    if stopPos==0:
        stopPos = len(seq)
    stopPos = min(stopPos,len(seq)-len(pattern))
    #print("seq:", seq)
    #print("pattern:", pattern)
    #print("startPos", startPos, "stopPos", stopPos)
    if startPos>=stopPos:
        print("Start position must be less than stop position.")
        sys.exit()

    if numMatches==-1:
        numMatches = len(seq)-len(pattern)
 
    counter = startPos
    countMatches = 0
    patternPositions = []
    while counter < stopPos and countMatches<numMatches:
        thispattern = seq[counter:counter+len(pattern)]
        if thispattern==pattern:
            patternPositions.append(counter+1) #print(counter+1) #returning 1-based position
            countMatches+=1
        counter+=1
        
    return patternPositions #returning an array of (1-based) positions where the pattern occurs in the given sequence

if __name__ == "__main__":
    #Script to tell me where a given sequence pattern occurs in a given larger sequence.
        #Derived from scanForORFs.py
    #It takes a sequence, by itself or in fasta format, and gives the relative position(s) of the pattern in the sequence.
    #Limitations: Input sequence has to all be on one line.
    #TO RUN: python scanForORFs.py seq.txt/fa
    #Can also provide a (1-based) start position to start scanning or an end position to stop scanning
    #Returns the (1-based_ relative position of each occurence of the pattern in the sequence.

    startPos = 0
    stopPos = 0
    inputName = ""
    pattern = ""
    numMatches = -1
    #############Begin command line options code
    argv = sys.argv[1:] #grabs all the arguments

    def usage():
        print('''python scanForSeq.py -i input_sequence -p sequence_pattern_to_search_for [-a, -b] > output
        Returns the (1-based) position of occurences of pattern -p in input sequence -i.
        Options:
        -h, --help

        REQUIRED
        -i, --input         The input sequence or sequence file.
        -p, --pattern       The sequence pattern to find in input.

        OPTIONAL
        -a, --start         Position to start scanning.
        -b, --stop			Position to stop scanning.
        -n, --num			Stop looking after n matches found.
                            If both -b and -n used, will stop at whichever occurs first.
        ''')

    if len(argv)==0:
        usage()
        sys.exit()

    try:
        opts, args = getopt.getopt(argv, "hi:p:a:b:n:", ["help","input=", "pattern=","start=","stop=", "num="])
    except getopt.GetoptError: #option help input log and sqrt
        usage()
        sys.exit(2)

    for opt, arg in opts:           
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            inputName = arg

        elif opt in ("-p", "--pattern"):
            pattern = arg.upper().replace('U','T')

        elif opt in ("-a", "--start"):
            startPos=int(arg)-1

        elif opt in ("-b", "--stop"):
            stopPos=int(arg)-1

        elif opt in ("-n", "--num"):
            numMatches=int(arg)

    if pattern=="" or inputName=="":
        usage()
        sys.exit()
    source = "".join(args)
    if len(source)>0:
        print("WARNING: Unused options", source)
    ################# END command line options code
    positions = scanForSeq(inputName,pattern,startPos,stopPos,numMatches)
    for i in positions:
        print(i)
