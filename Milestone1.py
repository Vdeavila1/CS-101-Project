{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "c8269f7b-b23e-4af1-9119-ec577cec1da1",
   "metadata": {},
   "outputs": [],
   "source": [
    "def reverse_complement(dna):\n",
    "## This function should accept a string representing a DNA sequence and return the reverse complement as a new string. The reverse complement of a DNA string is formed by reversing the entire string, and then taking the complement of each nucleotide.\n",
    "    dnaComp = ''\n",
    "    dnaList = list( dna.upper() )\n",
    "    dnaList.reverse()\n",
    "    for letter in dnaList:\n",
    "        if letter == 'A':\n",
    "            dnaComp += 'T'\n",
    "        elif letter == 'C':\n",
    "            dnaComp += 'G'\n",
    "        elif letter == 'G':\n",
    "            dnaComp += 'C'\n",
    "        elif letter == 'T':\n",
    "            dnaComp += 'A'\n",
    "    return dnaComp"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "c1bfcca1-78d6-47aa-b703-0a04594c340b",
   "metadata": {},
   "outputs": [],
   "source": [
    "def GC_content(dna_list):\n",
    "##This function should accept a list of DNA strings, and return the index of the DNA string with the highest GC-content and its GC-content percentage as a tuple.\n",
    "##The GC-content of a DNA string is the percentage of nucleotides in the string that are \"C\" or \"G\".\n",
    "    contentList = []\n",
    "    for dnaString in dna_list:\n",
    "        countGC = 0\n",
    "        countAT = 0\n",
    "        for i in dnaString.upper():\n",
    "            if letter == 'G' or letter == 'C':\n",
    "                countGC += 1\n",
    "            elif letter == 'A' or letter == 'T':\n",
    "                countAT += 1\n",
    "        contentList.append( countGC / ( countGC + countAT ) * 100 )\n",
    "    maximum = contentList[ 0 ]\n",
    "    for content in contentList:\n",
    "        if content > maximum:\n",
    "            maximum = content\n",
    "    return ( contentList.index( maximum ), maximum )\n",
    "                \n",
    "            "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "8219f3ea-a568-4e29-ae77-57523c22acc3",
   "metadata": {},
   "outputs": [],
   "source": [
    "def rna2codon(rna):\n",
    "## This function should accept a string representing an RNA sequence, and return the corresponding amino acid string, as transcribed by this codon table. \n",
    "## You do not need to transcribe the stop codon.\n",
    "    genetic_code = {\n",
    "        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',\n",
    "        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',\n",
    "\n",
    "        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',\n",
    "        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',\n",
    "\n",
    "        'UAU': 'Y', 'UAC': 'Y', 'UAA': 'STOP', 'UAG': 'STOP',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',\n",
    "        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',\n",
    "\n",
    "        'UGU': 'C', 'UGC': 'C', 'UGA': 'STOP', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',\n",
    "        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',\n",
    "    }\n",
    "    allowed_codons = set('ACGU')\n",
    "    if rna in genetic_code.keys():\n",
    "        return genetic_code[rna]\n",
    "    else: \n",
    "        return 'Invalid'"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "880e242a-acfa-46ca-ab81-f16b46a6592b",
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "b6286c18-7801-48b5-866a-934113d2e0ce",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.8"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
