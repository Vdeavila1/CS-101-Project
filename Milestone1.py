def reverse_complement(dna):
## This function should accept a string representing a DNA sequence and return the reverse complement as a new string. The reverse complement of a DNA string is formed by reversing the entire string, and then taking the complement of each nucleotide.\n",
    dnaComp = ''
    dnaList = list( dna.upper() )
    dnaList.reverse()
    for letter in dnaList:
        if letter == 'A':
            dnaComp += 'T'
        elif letter == 'C':
            dnaComp += 'G'
        elif letter == 'G':
            dnaComp += 'C'
        elif letter == 'T':
            dnaComp += 'A'
    return dnaComp

def GC_content(dna_list):
##This function should accept a list of DNA strings, and return the index of the DNA string with the highest GC-content and its GC-content percentage as a tuple.\n",
##The GC-content of a DNA string is the percentage of nucleotides in the string that are \"C\" or \"G\".\n",
    contentList = []
    for dnaString in dna_list:
        countGC = 0
        countAT = 0
        dnaString = dnaString.upper()
        for letter in dnaString:
            if letter == 'G' or letter == 'C':
                countGC += 1
            elif letter == 'A' or letter == 'T':
                countAT += 1
        contentList.append( countGC / ( countAT + countGC ) * 100 )
    maximum = contentList[ 0 ]
    for content in contentList:
        if content > maximum:
            maximum = content
    return ( contentList.index( maximum ), maximum )
  

def rna2codon(rna):
    geneticCode = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
 
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
 
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '', 'UAG': '',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
 
        'UGU': 'C', 'UGC': 'C', 'UGA': '', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    output = ''
    threes = [rna[x:x+3] for x in range(0, len(rna), 3)]  # Splits rna into groups of 3
    for chunk in threes:                                                           # adds each rna protein to final list
        if geneticCode[chunk] == '':                                     # stops if it hits a stop codon
            break
        output += geneticCode[chunk]
    return output
rna = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
print(rna2codon(rna))
 
def s(dna):
    final_dict = {}
    total_ncltde_count = []
    nucleotides = ['A','C','G','T']
    dna = dna.upper()
    for c in nucleotides: #counts how many instances of each nucleotide
        ncltde_count = dna.count(c)
        total_ncltde_count += [ncltde_count]
    for num in range(len(nucleotides)): # adds the number to dict with its respective nucleotide
        final_dict[nucleotides[num]] = total_ncltde_count[num]
    return(final_dict)


def mendels_law(hom, het, rec):
    pop = hom + het + rec # total population
    d = hom
    b = het
    r = rec
    prob = (8*r*d+8*d*b+4*r*b+3*(b**2-b)+4*(d**2-d))/(4*(pop**2-pop))
    return prob

def hamming_dist(dna1,dna2):
    count = 0
    for i in range(len(dna1)): #for each nucleotide in dna1
        if dna1[i] is dna2[i]: #check if equal to dna2's corresponding nucleotide 
            continue           #if so, move on
        else:
            count += 1         #if not equal, add to counter for total # of differences
    return count


def count_dom_phenotype(genotypes):
    count = 0
    for n in genotypes[:3]:      # takes the first 3 ints from the list
        count += n               # add to count (Anytime there is AA parent, all childs have 1 dom allele)
    count += genotypes[3] * 0.75 # 4th int is Aa * Aa which is 0.75 dom alleles
    count += genotypes[4] * 0.5  # 5th int is Aa * aa whcih is 0.5 alleles
    return (count*2)             # mult by 2 since 2 children per pair

def splice_rna(dna, intron_list):
##This function should accept a string representing a DNA sequence and a list of strings representing introns. 
##The process of transcribing DNA into RNA involves translating the DNA to RNA and then performing RNA splicing, where the sequence is chopped into smaller segments called introns and exons
##Introns are segments of the gene not used for protein translation, so they should be removed from the sequence.
##Exons are the remaining segments, which are then transcribed sequentially into a protein string. splice_rna() should return a protein string that results from transcribing and translating the exons of the given string.
    for x in intron_list:
        dna = dna.replace( x, '' )
    rna = dna2rna( dna )    ##Converts dna inputs to rna using dna2rna function
    codon = rna2codon( rna )    ##Converts rna input to codons using rna2codon function
    return codon

def dna2rna(dna):
##This function should accept a string representing a DNA sequence and return the transcribed RNA string. 
##DNA sequences are transcribed to RNA sequences by replacing all occurrences of the "T" nucleotides with "U"
    rna = ''
    for symbol in dna:
        if symbol == 'T':
            rna = rna + 'U'
        elif symbol == 'A':
            rna = rna + 'A'
        elif symbol == 'G':
            rna = rna + 'G'
        else:
            rna = rna + 'C'
    return rna

def fibonacci_rabbits(n,k):
##This function should accept two integers n and k and calculate the total number of rabbit pairs that will be present after nn months, given that every pair of mating rabbits produces k rabbit pairs in their next litter (in the next month).
    month1=1
    month2=1
    for i in range(n-1):
        month2, month1 =month1, month1+(month2*k)
    return month2
   
def source_rna(protein):
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    i = 0
    total = 1
    for i in range(len(protein)):#goes through each letter of the protein string
        occur = 0
        for key in genetic_code.keys():#keys of genetic_code
            if genetic_code[key] == protein[i]: #If the keys in genetic_code = the value of the ith place in protein then adds 1 to occur
                occur+=1
        total = total * occur
        i+=1
    return total*3
print(source_rna("CDMA"))

def locate_substring(dna_snippet,dna):
    occur_list = []
    count = 0
    i=0
    while (i<len(dna)):#loops to find where each occurence is
        count = dna.find(dna_snippet,count)#finds first occurence of dna_snippet and stores in count
        if count == -1:
            break
        else:
            occur_list.append(count)
        i+=1
        count+=1
    return occur_list
dnas = "ATAT"
dna0 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA"
print(locate_substring(dnas,dna0))
