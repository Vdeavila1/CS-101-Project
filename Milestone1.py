
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
            for i in dnaString.upper():
                if letter == 'G' or letter == 'C':
                    countGC += 1
                elif letter == 'A' or letter == 'T':
                    countAT += 1
            contentList.append( countGC / ( countGC + countAT ) * 100 )
        maximum = contentList[ 0 ]
        for content in contentList:
            if content > maximum:
                maximum = content
        return ( contentList.index( maximum ), maximum )
  
    def rna2codon(rna):\n",
    ## This function should accept a string representing an RNA sequence, and return the corresponding amino acid string, as transcribed by this codon table. \n",
    ## You do not need to transcribe the stop codon.\n",
        genetic_code = {\n",
            'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',\n",
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',\n",
    
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',\n",
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',\n",
                        
            'UAU': 'Y', 'UAC': 'Y', 'UAA': 'STOP', 'UAG': 'STOP',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',\n",
            'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',\n",
    
            'UGU': 'C', 'UGC': 'C', 'UGA': 'STOP', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',\n",
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',\n",
        
        allowed_codons = set('ACGU')\n",
        if rna in genetic_code.keys():\n",
            return genetic_code[rna]\n",
        else: \n",
            return 'Invalid'"
   
  
  
   
