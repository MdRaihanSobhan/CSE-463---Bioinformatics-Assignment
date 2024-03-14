import random
from datetime import datetime

def symbolToNumber(symbol): # symbol is a string of length 1 
	if symbol == "A":
		return 0
	if symbol == "C":
		return 1
	if symbol == "G":
		return 2
	if symbol == "T":
		return 3

def numberToSymbol(x): # x is an integer 
	if x == 0:
		return "A"
	if x == 1:
		return "C"
	if x == 2:
		return "G"
	if x == 3:
		return "T"

def profileRandom(k, profile, text): # k is an integer, profile is a list of lists of floats, text is a string 
    probs = [] 
    for i in range(0,len(text) - k +1): # for each kmer in text 
        prob = 1.0
        pattern = text[i:i+k] # pattern is the kmer 
        for j in range(k): # for each symbol in the kmer
            l = symbolToNumber(pattern[j]) # l is the index of the symbol in the profile 
            prob *= profile[l][j] # prob is the product of the probabilities of the symbols in the kmer 
        probs.append(prob) # add the probability of the kmer to the list of probabilities
    r = myRandom(probs) # r is the index of the kmer with the highest probability 
    return r

def profileForm(motifs): # motifs is a list of strings
	k = len(motifs[0]) # k is the length of the first string in motifs
	profile = [[1 for i in range(k)] for j in range(4)] # initialize profile to a list of lists of length k  with all elements initialized to 1 
	for x in motifs: # for each string in motifs
		for i in range(len(x)): # for each symbol in the string
			j = symbolToNumber(x[i]) # j is the index of the symbol in the profile 
			profile[j][i] += 1 # add 1 to the count of the symbol in the profile
	for x in profile: # for each list in profile
		for i in range(len(x)): # for each element in the list
			x[i] = x[i]/len(motifs) # divide the count of the symbol by the number of strings in motifs
	return profile

def consensus(profile): # profile is a list of lists of floats
	str = ""
	for i in range(len(profile[0])): # for each element in the first list in profile
		max = 0 
		loc = 0
		for j in range(4): # for each element in the profile
			if profile[j][i] > max: # if the element in the profile is greater than the current max
				loc = j # set the location to the index of the element in the profile
				max = profile[j][i] # set the max to the element in the profile
		str+=numberToSymbol(loc) # add the symbol corresponding to the location to the string 
	return str

def score(motifs): # motifs is a list of strings
	profile = profileForm(motifs) # profile is a list of lists of floats 
	cons = consensus(profile) # cons is a string 
	score = 0 
	for x in motifs: # for each string in motifs
		for i in range(len(x)): # for each symbol in the string
			if cons[i] != x[i]: # if the symbol in the consensus is different from the symbol in the string
				score += 1 # add 1 to the score
	return score

def myRandom(dist): # dist is a list of floats
    s = 0.0
    for x in dist: # for each element in dist
        s+= x # add the element to s
    i = random.random() # i is a random number between 0 and s
    partial = 0.0
    for x in range(len(dist)): # for each element in dist
        partial += dist[x] # add the element to partial
        if partial/s >= i: # if the partial sum divided by s is greater than or equal to i
            return x

def gibbsSampler(dna, k, t, n): # dna is a list of strings, k,t and n are integers
    bestMotifs = []
    motifs = []
    for x in range(t): # for each string in dna
        i = random.randint(0, len(dna[x])-k) # i is a random integer between 0 and the length of the string minus k
        motifs.append(dna[x][i:i+k]) # add the k-mer starting at index i to the motifs list
    bestMotifs = motifs[:] # bestMotifs is a copy of motifs
    for i in range(n): # for each iteration
        j = random.randint(0,t-1) # j is a random integer between 0 and t-1
        profile = profileForm(motifs[:j] + motifs[j+1:]) # profile is a list of lists of floats
        r = profileRandom(k, profile, dna[j]) # r is a random k-mer from the profile
        motifs[j] = dna[j][r:r+k] # set the jth element of motifs to the k-mer starting at index r
        if score(motifs) < score(bestMotifs): # if the score of motifs is less than the score of bestMotifs
            bestMotifs = motifs[:] # set bestMotifs to a copy of motifs
    return bestMotifs

start = datetime.now()
with open("yst08r.txt") as file:
    #k=Motif Length # We will run for 8 to 15
    #t=  indicating how many times to read the genetic algorithm.
    #N Interation count
    k, t, n = [int(x) for x in file.readline().split()] 
    #print(k,t,n)
    dna = []
    for line in file:
        dna.append(line.rstrip()) 

best = gibbsSampler(dna, k, t, n) 
s = score(best)

for x in range(20): # for each iteration of the gibbs sampler 
    sample = gibbsSampler(dna, k, t, n)
    if score(sample) < s: # if the score of the sample is less than the score of the best
        s = score(sample) # set s to the score of the sample
        best = sample[:] # set best to a copy of the sample
for b in best: 
	print(b)


end = datetime.now()
print("Runtime : ",end,start)