"""
Based on Material taught from Bioinformatics Algorithms by Phillip Compeau, some functions are based on pseudocode found in the text
All code is written by Michael Taylor
OriC Finder, started Feb 16, 2022
"""


def substring(text, i, k):
    return text[i:i+k]


def findMinSkew(genome):
    """ Returns position of minimum skew, skew = # of G's - # of C's

    Args:
        genome (_type_): DNA sequence as string. e.g., "ACTGATGC"

    Returns:
        list: list of positions/indices in genome of lowest skew
    """
    minSkewIndices = []
    curSkew = 0  # skew at current position in genome
    minSkew = 0  # current minimum skew
    numGs = 0
    numCs = 0

    # letter represents the index of the letter in the genome in this loop
    for genomePosition in range(len(genome)):
        # Keep running total of # of G & C occurences
        if genome[genomePosition].upper() == "G":
            numGs += 1
        if genome[genomePosition].upper() == "C":
            numCs += 1
        # Update skew at current position
        curSkew = numGs - numCs
        # Record absolute minimums
        if curSkew < minSkew:
            minSkew = curSkew
            minSkewIndices = []
            minSkewIndices.append(genomePosition + 1)
        if curSkew == minSkew:
            for indice in minSkewIndices:
                if genomePosition-indice > 10:  # makes sure nearby indices arent included
                    minSkewIndices.append(genomePosition + 1)
    # Return list of indices
    return minSkewIndices


def findOri(genome):
    # Step 1: Find position of #G-#C skew's absolute minimum
    # list of minimum skew indices/position(s) in the genome
    minSkew = findMinSkew(genome)

    return "Ori Found?"


ecoliGenomePath = r"c:\Users\micha\Desktop\OriC Finder\E_coli.txt"
ecoliGenomeFile = open(ecoliGenomePath)
ecoliGenome = ecoliGenomeFile.read()

findOri(ecoliGenome)
