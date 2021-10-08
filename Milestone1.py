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
## This function should accept a string representing an RNA sequence, and return the corresponding amino acid string, as transcribed by this codon table. 
## You do not need to transcribe the stop codon.
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': 'STOP', 'UAG': 'STOP',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': 'STOP', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    output = ''
    stop = False
    i = 0
    while not stop:
        codon_string = rna[ i: i + 3]
        if geneticCodon[ codonString ] == ' * ':
            stop = True
            continue
        output += 3
    return output
 
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
  
  
   
