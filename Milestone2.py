from math import factorial
def find_splice( dna_motif, dna ):
  ##This function should accept as argument two strings representing a DNA motif dna_motif (a subsequence of dna) and a DNA sequence dna,
  ##and return the index positions of the motif characters for the first occurrence of the substring within dna.
  ##A subsequence is a collection of symbols contained in a specific order but not necessarily contiguously.
  ##This means that the pattern may be interrupted by other characters. 
  ##(HINT: recall that Python is 0-indexed, meaning we start counting positions at 0...)
    character = dna_motif[ 0 ]
    index = 0
    positions = []
    for x in range( len( dna ) ):
        if (character == dna[ x ]):
            positions.append( x )
            index = index + 1
            if index >= len( dna_motif ):
                break
            character = dna_motif[ index ]
    if index != len( dna_motif ):
        return []
    return positions

  

def assemble_genome( dna_list ):
    ## This function should accept a list of DNA strings and return the shortest superstring containing all given DNA strings.
    ## A superstring is a string that contains each of the smaller provided strings as a substring. 
    ## Recall, a substring is a contiguous sequence of characters within a string. Assume that the strings will all be the same length.
    dictionary = {}
    for i in range( len( dna_list ) ):
        for j in range( len( dna_list ) ):
            if i != j:
                a = 0
                for k in range( 1, min( len( dna_list [i] ), len( dna_list [j] ) ) ):
                    if dna_list [j] [:k] == dna_list [i] [-k:]:
                        a = k
                dictionary [ ( i, j ) ] = a
    if max( dictionary.values() ) == 0:
        return ''.join( dna_list )
    else:
        answer = ''.join( dna_list )
        length = len( answer )
        return_list = []
        for i, word in enumerate( dna_list ):
            pack = set( range( len( dna_list ) ) )
            pack.remove( i )
            return_list.append( ( word, i, pack) )
        while return_list:
            v1, pl, left = return_list.pop()
            if len( v1 ) < length:
                if not left:
                    answer = v1
                    length = len( answer )
                else:
                    pack = [ [ dictionary[ pl, index ], index] for index in left ]
                    pack.sort()
                    for a, index in pack:
                        v2 = v1 + dna_list[ index ][ a: ]
                        last = set( left )
                        last.remove( index )
                        return_list.append( ( v2, index, last ) )
        return answer
      
def shared_motif(dna_list):
    String=dna_list[0]
    StringLength=len(String)
    Substrings=[]
    Result = ""
    for i in range(StringLength):
        for n in range(i+1,StringLength+1):
            #different lengths of substrings
            test=String[i:n]
            IsIn=True
            for m in range(1,len(dna_list)):
                #looks through all the other strings in the list
                if test not in dna_list[m]:
                    IsIn=False
                    break #break out of loop if substring is not in string
            if IsIn==True:
                Substrings.append(test)
    for Sub in Substrings:
        if len(Sub)>len(Result): #Finds longest substring
            Result=Sub
    return Result

print(shared_motif(["GATTACA","TAGACCA","ATACA"]))
print(shared_motif(["ATATACA", "ATACAGA", "GGTATACA"]))


# In[59]:


import math
def perfect_match(rna):
    rna = rna.upper()
    Acount = rna.count('A')
    Ucount = rna.count('U')
    Gcount = rna.count('G')
    Ccount = rna.count('C')
    AUcount = 0
    GCcount = 0
    count = 0
    if Acount == Ucount:
        AUcount = math.factorial(Acount)
    if Gcount == Ccount:
        GCcount = math.factorial(Gcount)
    count = GCcount*AUcount
    return count
testrna = "AAACCCGGGUUU"
print(perfect_match(testrna))


def reverse_complement(dna):
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



def rev_palindrome(dna):
    ans = []
    
    for x in range(len(dna)):           #iterate over each letter in the string (x = index)
    
        for i in range(4,12):   #Loop using slice range from 4 to 12 characters
            if (x + i) > len(dna):   #skip if index range oversteps bounds
                continue
            
            chunk = dna[ x : x + i]       #set chunk using each length 4 to 12, starting at current index
            revComp = reverse_complement( chunk )    #find reverse compliment of current chunk

            if revComp == chunk:        #if reverse compliment equals original chunk
                ans.append((x,i))    #add to list the position and length

    return ans
  
  import math

def random_genome(dna, gc_content):
    ans = []
    
    at_count = dna.count('A') + dna.count('T')
    gc_count = dna.count('G') + dna.count('C')

    for x in gc_content:
        gc_frequency = x/2
        at_frequency = (1-x)/2
    
        gc_result = ((gc_frequency)**gc_count)
        at_result = ((at_frequency)**at_count)
        
        ans.append( math.log10(gc_result * at_result) )
    
    return ans


def get_edges(dna_dict):
    dna_strings = list(dna_dict.values())   # list of all DNA strings
    prefixes = []
    suffixes = []
    ans = []
    
    rev_dict = { v:k for k,v in dna_dict.items() }  #reverses inputted dictionary so we can change strings to identifiers

    
    for e in dna_strings:               #makes a list of suffixes & prefixes
        prefixes.append( e[0:3] )       #first 3 letters
        suffixes.append( e[-1:-4:-1] )  #last 3 letters
    
    for i in range(len(dna_strings)):   #iterate through each dna string

        suf = suffixes[i]                 #set suffix for current string
        
        for x in range(len(suffixes)):                  #Iterate for each suffix
            if i == x:                                  #skip when the prefix and suffix come from the same string
                continue
            elif suf == prefixes[x]:                      #find suffixes that match prefixes
                preID = rev_dict.get(dna_strings[i])    #convert string with suffix to identifier
                suffID = rev_dict.get(dna_strings[x])   #convert string with prefix to identifier
                ans.append((preID,suffID))        # add IDs as a tuple to final list
    
    return ans
