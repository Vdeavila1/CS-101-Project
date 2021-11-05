def find_splice( dna_motif, dna ):
  ##This function should accept as argument two strings representing a DNA motif dna_motif (a subsequence of dna) and a DNA sequence dna,
  ##and return the index positions of the motif characters for the first occurrence of the substring within dna.
  ##A subsequence is a collection of symbols contained in a specific order but not necessarily contiguously.
  ##This means that the pattern may be interrupted by other characters. 
  ##(HINT: recall that Python is 0-indexed, meaning we start counting positions at 0...)
    position = []
    a = -1
    for i in dna_motif:
        x = dna.find( i, a + 1 )
        position.append( a )
    if -1 in position:
        return []
    return position
  
from math import factorial
def assemble_genome( dnaList ):
    ## This function should accept a list of DNA strings and return the shortest superstring containing all given DNA strings.
    ## A superstring is a string that contains each of the smaller provided strings as a substring. 
    ## Recall, a substring is a contiguous sequence of characters within a string. Assume that the strings will all be the same length.
    dictionary = {}
    for i in range( len( dnaList ) ):
        for j in range( len( dnaList ) ):
        if i != j:
            x = 0
            for k in range( 1, min( len( dnaList [i] ), len( dnaList [j] ) ) ):
                if dnaList [j] [:k] == dnaList [i] [-k:]:
                    x = k
                dictionary [ ( i, j ) ] = x
    if max( dictionary.values() ) == 0:
        return ''.join( dnaList )
    else:
        answer = ''.join( dnaList )
        length = len( answer )
        return_list = []
        for i,word in enumerate( dnaList ):
            pack = set( range( len( dnaList ) ) )
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
                    for n, index in pack:
                        v2 = v1 + dnaList[ index ][ n: ]
                        last = set( left )
                        last.remove( index )
                        return_list.append( ( v2, index, last ) )
        return answer
