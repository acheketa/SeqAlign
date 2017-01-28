#! /usr/bin/env python
'''
Pairwise sequence alignment using BLOSUM62

Niema Moshiri 2017
'''
from common import BLOSUM62 as M
from common import indelToGap
from common import parseFASTA
import argparse,time

# perform pairwise alignment
def pairAlign(seqs, gap):
    if len(seqs) != 2:
        print("ERROR: pairAlign must take in exactly 2 sequences")
        exit(-1)
    sID,tID = sorted(seqs.keys())
    s,t = seqs[sID],seqs[tID]
    S = [[0 for j in range(len(t)+1)] for i in range(len(s)+1)]
    B = [[None for j in range(len(t)+1)] for i in range(len(s)+1)]
    for i in range(1, len(s)+1):
        S[i][0] = i*gap
        B[i][0] = 'i'
    for j in range(1, len(t)+1):
        S[0][j] = j*gap
        B[0][j] = 'j'
    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            options = {'ij': S[i-1][j-1] + M[s[i-1]][t[j-1]], 'i': S[i-1][j] + gap, 'j': S[i][j-1] + gap}
            S[i][j] = max(options.values())
            for o in options:
                if options[o] == S[i][j]:
                    B[i][j] = o
                    break
    i,j = len(s),len(t)
    out = {sID:'', tID:''}
    while i+j != 0:
        bt = B[i][j]
        if 'i' in bt:
            out[sID] += s[i-1]
            i -= 1
        else:
            out[sID] += '-'
        if 'j' in bt:
            out[tID] += t[j-1]
            j -= 1
        else:
            out[tID] += '-'
    out[sID] = out[sID][::-1]
    out[tID] = out[tID][::-1]
    return S[len(s)][len(t)], out

# main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',  required=True,  type=argparse.FileType('r'), help="Input FASTA File")
    parser.add_argument('-o', '--output', required=True,  type=argparse.FileType('w'), help="Output FASTA File")
    parser.add_argument('-p', '--indel',  required=False, type=float, default=0.1,     help="Indel Ratio")
    args = parser.parse_args()
    start_time = time.time()
    print("User-specified indel rate: " + str(args.indel))
    gap = indelToGap(args.indel)
    print("Indel rate converted to gap rate: " + str(gap))
    print("Parsing input FASTA file...")
    seqs = parseFASTA(args.input)
    print("Computing 2D dynamic programming matrix...")
    score,aln = pairAlign(seqs, gap)
    print("Writing output alignment...")
    for ID in sorted(aln.keys()):
        args.output.write('>')
        args.output.write(ID)
        args.output.write('\n')
        args.output.write(aln[ID])
        args.output.write('\n')
    print("--- Completed in %s seconds ---" % (time.time() - start_time))