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
        if (character == dna[ i ]):
            positions.append( x )
            index = index + 1
            if index >= len( dna_motif ):
                break
            character = dna_motif[ index ]
    if index != len( dna_motif ):
        return []
    return positions

  
from math import factorial
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
