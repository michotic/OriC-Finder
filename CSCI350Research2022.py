# These functions were initially made for coding tests while studying Ch.1 of Bioinformatics Algorithms in CSCI350 at St.FX
# In this file I have modified them to make necessary changes for finding better ori candidates.
# E. coli will be the initial DNA sequence used to optimize my "ori finding" algorithm
#
#
# 1)Slight variations of the target k-mer of <= d hammingDistance
#
from copy import copy
from operator import ge


def substring(text, i, k):
    return text[i:i+k]


# Returns int:count of occurences of str:pattern in str:genome
#Complexity: O(n)
def patternCount(genome, pattern):
    count = 0
    for i in range(len(genome) - len(pattern)):
        if substring(genome, i, len(pattern)) == pattern:
            count = count + 1
    return count

# Returns list of most frequent patterns of int:k length in str:genome
# Complexity: O(n^2)


def freqWords(genome, k):
    freqPatts = []
    counts = []
    for i in range(len(genome) - k):  # O(n)
        pattern = genome[i:i+k]
        # O(n) function called O(n) times -> O(n^2)
        counts.append(patternCount(genome, pattern))

    maxCount = max(counts)

    for i in range(len(genome) - k):  # O(n)
        if counts[i] == maxCount:
            freqPatts.append(genome[i:i+k])
    return list(set(freqPatts))

# Returns dict of most frequent patterns of int:k length in str:genome


# Prints list of most frequent k-mers through lengths minK to maxK


def makeKmerTable(genome, minK, maxK):
    for i in range(minK, maxK + 1):
        freqTable = freqWordsPlus(genome, i)
        print(" ")
        print("frequent k-mers of length k = " +
              str(i) + ", count = " + str(len(freqTable)))
        print(freqTable)
        print(" ")
        print("_______________________________________________________________")

# outputs the reverse complementary strand of the given DNA string


def reverseDNA(genome):
    dnaPairings = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }
    backwardsDNA = genome[::-1]
    complementaryStrand = ""
    for nucleotide in backwardsDNA.upper():
        complementaryStrand += dnaPairings.get(nucleotide)
    return complementaryStrand


# Finds clumps of k-mers that occur >= t times within range L of each other within genome
def findClumps(genome, k, L, t):
    patterns = []
    n = len(genome)
    for i in range(n-L):
        window = genome[i:i+L]
        freqMap = freqTable(window, k)
        for kmer in freqMap.keys():
            if freqMap[kmer] >= t:
                patterns.append(kmer)
    patterns = list(set(patterns))
    return patterns

# Calculates skew of G to C in genome (# of G's - # of C's)


def skew(genome):
    curSkew = 0
    for letter in genome:
        if letter.upper() == "G":
            curSkew += 1
        if letter.upper() == "C":
            curSkew -= 1
    return curSkew

# returns a list of all indices in the genome that are a point of lowest skew


def findMinSkew(genome):
    minSkewIndices = []
    curSkew = 0
    minSkew = 0

    # letter represents the index of the letter in the genome in this loop
    for letter in range(len(genome)):
        if genome[letter].upper() == "G":
            curSkew += 1
        if genome[letter].upper() == "C":
            curSkew -= 1
        if curSkew < minSkew:
            minSkew = curSkew
            minSkewIndices = []
            minSkewIndices.append(letter + 1)
        if curSkew == minSkew:
            minSkewIndices.append(letter + 1)
    return list(set(minSkewIndices))


def hammingDistance(kmer1, kmer2):
    hamDist = 0
    for i in range(min([len(kmer1), len(kmer2)])):
        if kmer1[i] != kmer2[i]:
            hamDist += 1
    return hamDist


# Returns a collection of space-separated integers specifying all starting positions where Pattern appears as a substring of Genome
def findPatternPositions(pattern, genome):
    indices = []
    for i in range(len(genome) - len(pattern)):
        if genome[i:i+len(pattern)] == pattern:
            indices.append(i)
    return indices


def findSimilarPatternPositions(pattern, genome, d):
    positions = []  # list of positions for similar patterns in genome sequence
    k = len(pattern)

    for i in range(len(genome) - len(pattern) + 1):
        kmer = substring(genome, i, k)
        similar = False

        # Check 1: Near match?
        if hammingDistance(pattern, kmer) <= d:
            similar = True
        # Check 2: Is it a near match for the reverse complement?
        if hammingDistance(reverseDNA(pattern), reverseDNA(kmer)) <= d:
            similar = True

        if similar:
            positions.append(i)
    return positions


'''
IterativeNeighbors(Pattern, d)
        Neighborhood ← set consisting of single string Pattern
        for j = 1 to d
            for each string Pattern’ in Neighborhood
                add ImmediateNeighbors(Pattern') to Neighborhood
                remove duplicates from Neighborhood
        return Neighborhood
'''


def immediateVariations(pattern):
    neighborhood = [pattern]
    for i in range(len(pattern)):
        symbol = pattern[i]
        otherSymbols = ["A", "T", "G", "C"]
        if symbol in otherSymbols:
            otherSymbols.remove(symbol)
        for newSymbol in otherSymbols:
            newVariation = list(pattern)
            newVariation[i] = newSymbol
            neighborhood.append("".join(newVariation))

    return neighborhood


def similarPatternsList(pattern, d):
    # Returns list of all possible similar patterns with at most d differences
    neighborhood = [pattern]
    n = len(neighborhood)
    for j in range(d):
        for kmer in range(n):
            neighbors = immediateVariations(neighborhood[kmer])
            for thing in neighbors:
                if thing not in neighborhood:
                    neighborhood.append(thing)
            # print(neighbors)
        n = len(neighborhood)
    return neighborhood


def countNearMatches(genome, pattern, d):
    potentialMatches = similarPatternsList(pattern, d)
    k = len(pattern)
    count = 0
    for position in range(len(genome) - k + 1):
        if substring(genome, position, k) in potentialMatches:
            count += 1
    return count


def similarPatternCount(genome, pattern, d):
    count = 0
    k = len(pattern)
    for position in range(len(genome) - k + 1):
        kmer = substring(genome, position, k)
        if hammingDistance(pattern, kmer) <= d:
            count += 1
    return count


def freqSimilarWords(genome, k, d):
    freqPatts = []
    maxFreq = 0

    for i in range(len(genome)-k+1):
        kmer = substring(genome, i, k)
        count = similarPatternCount(genome, kmer, d)
        firstOccurence = ""

        if count >= maxFreq:
            newKmer = True

            if count > maxFreq:
                freqPatts.clear()

            for patt in freqPatts:  # <<<?? this doesnt have anythig yet duhgoejsfiufbnwijf
                if hammingDistance(kmer, patt) <= d:
                    newKmer = False
                    firstOccurence = patt

            if newKmer:
                freqPatts.append(kmer)
                maxFreq = count
            else:
                freqPatts.append(firstOccurence)

    return list(set(freqPatts))


print(freqSimilarWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1))

ecoliGenomePath = r"c:\Users\micha\Desktop\OriC Finder\E_coli.txt"
ecoliGenomeFile = open(ecoliGenomePath)
ecoliGenome = ecoliGenomeFile.read()
