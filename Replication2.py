# Copy your PatternCount function from the previous step below this line

#Text = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"

#Pattern = "TGATCA"

#Text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"

#k = 4

#TargetURL=http://bioinformaticsalgorithms.com/data/realdatasets/Replication/Vibrio_cholerae.txt
#import urllib.request # the lib that handles the url stuff

#with urllib.request.urlopen("http://bioinformaticsalgorithms.com/data/realdatasets/Replication/Vibrio_cholerae.txt") as response:
    #Text = response.read()

#Text = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"



k = 3

#Pattern= "CTTGATCAT"
Pattern = "AAA"
Text = "TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT"

ecoli="GATATATGCATATACTT"
symbol = "C"

Genome = "CATGGGCATCGGCCATACGCC"
Genome2 = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"


def FrequentWords(Text, k):
    FrequentPatterns = [] # output variable
    # your code here
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternsNoDuplicates = remove_duplicates(FrequentPatterns)
    return FrequentPatternsNoDuplicates

def remove_duplicates(Items):
    ItemsNoDuplicates = [] # output variable
    # your code here
    seen = set()
    for x in Items:
        if x not in seen:
            ItemsNoDuplicates.append(x)
            seen.add(x)
    return ItemsNoDuplicates




def CountDict(Text, k):
    Count = {} 
    # fill in the rest of the CountDict function here.
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Pattern, Text)
    return Count


def FasterSymbolArray(Genome, symbol):
    array = {}
    # your code here
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array
    
#Using a list
def Skew1(Genome):
    skew = [0] #initializing the list
    # your code here
    n = len(Genome)
    for i in range(0, n):
        if Genome[i] == "G":
            skew.append(skew[i]+1)
        elif Genome[i] == "C":
            skew.append(skew[i]-1)
        else:
            skew.append(skew[i])
    return skew


# Using a dictionary
def Skew(Genome):
    skew = {} #initializing the dictionary
    # your code here
    skew[0] = 0
    n = len(Genome)
    for i in range(1, n+1):
        if Genome[i-1] == "G":
            skew[i]= skew[i-1]+1
        elif Genome[i-1] == "C":
            skew[i]= skew[i-1]-1
        else:
            skew[i]= skew[i-1]
    return skew

def minimum_skew(Genome2):
    minimal_skew_ints = []


    skew = Skew(Genome2)
    m = min(skew.values())

    for i in range(0, len(Genome2)):
        current_skew = skew[i]
        if current_skew == m:
            minimal_skew_ints.append(i)
        elif current_skew < m:
            mimimal_skew = current_skew
            minimal_skew_ints = []
            minimal_skew_ints.append(i)

    return minimal_skew_ints
    

def HammingDistance(p, q):
    # your code here
    num= 0
    for i in range(0, len(p)):
        if p[i] != q[i]:
            num += 1
    return num
    
        


def PatternCount(symbol, Genome):
    count = 0
    for i in range(len(Genome)-len(symbol)+1):
        if Genome[i:i+len(symbol)] == symbol:
            count = count+1
    return count


def PatternMatch(Pattern, Genome):
    positions = []
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions



def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    # your code here
    for i in range(len(Genome) - len(Pattern) +1): #Loop over alignment
        match= True
        for j in range(len(Pattern)):  #Loop over characters
            if Genome[i+j] != Pattern[j]:
                match= False
                break
        if match:
           positions.append(i)
    return positions


def ApproximatePatternMatching(Pattern, Text, d):
    positions = [] # initializing list of positions
    # your code here
    for i in range(len(Text) - len(Pattern) + 1): # loop over alignment
        mismatch = 0
        match= True
        for j in range(len(Pattern)):  #loop over characters
            if Text[i+j] != Pattern[j]:
                #mismatch = HammingDistance(Pattern, Text)
                mismatch +=1
                if mismatch > d:
                    break
        if mismatch <= d:
           positions.append(i)
    return positions

## Not solve yet
def ApproximatePatternMatching2(Pattern, Text, d):
    positions = [] # initializing list of positions
    # your code here
    for i in range(len(Text) - len(Pattern) + 1): # loop over alignment
        mismatch = 0
        match= True
        #for j in range(len(Pattern)):  #loop over characters
            
        mismatch = HammingDistance(Pattern, Text)
                
        if mismatch > d:
            break
        if mismatch <= d:
           positions.append(i)
    return positions




def complement(Nucleotide):
    comp = '' #output variable
    # your code here
    complement = {'A':'T', 'C': 'G', 'G': 'C', 'T':'A'}
    bases = list(Nucleotide)
    bases =[complement[base] for base in bases]
    return ''.join(bases)


def ReverseComplement(Pattern):
    revComp = '' # output variable
    # your code here

    revComp = complement(Pattern[::-1])
    return revComp
            
   
 #   Genome = "GATATATGCATATACTT"
   # Pattern= "ATAT"

positions = PatternMatch(Pattern, Text)

Pattern1= "ATTCTGGA"
Text1 = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
d =3
pos_match= ApproximatePatternMatching(Pattern1, Text1, d)
pos_match2= ApproximatePatternMatching2(Pattern1, Text1, d)
skew1 = Skew(Genome)
skew2 = Skew1(Genome)
    
min_skew = minimum_skew(Genome2)

match = PatternCount(Pattern, Text)
print("the number of match is:\n")
print(match)

freq = FrequentWords(Text, k)
print("The most frequenct 10-mer is:\n")
print(freq)
print("The positions where the match is:\n")
print(positions)
print("Approximate match Pattern:\n")
print(pos_match)
print("Approximate match Pattern2:\n")
print(pos_match2)
print("FastSYmbolFunction:\n")
print(FasterSymbolArray(ecoli, symbol))
print("Skew genome:\n")
print(skew1)
print("Skew genome2:\n")
print(skew2)
print("Min Skew:\n")
print(min_skew)



# Now, set Text equal to the ori of Vibrio cholerae and Pattern equal to "TGATCA"


# Finally, print the result of calling PatternCount on Text and Pattern.
# Don't forget to use the notation print() with parentheses included!
