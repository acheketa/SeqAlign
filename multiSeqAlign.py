#! /usr/bin/env python
'''
Multiple sequence alignment using BLOSUM62

Niema Moshiri 2017
'''
from common import BLOSUM62 as M
from common import indelToGap
from common import parseFASTA
from pairSeqAlign import pairAlign
from tripleSeqAlign import tripleAlign
import argparse,time
try:
    import Queue as Q  # ver. < 3.0
except ImportError:
    import queue as Q

# merge single sequence with alignment using DP (which column does each symbol best fit)
def mergeAlign(aln, ID, seq, gap):
    keys = sorted(aln.keys())
    n = len(aln[keys[0]])
    m = len(seq)
    S = [[0 for j in range(m+1)] for i in range(n+1)]
    B = [[None for j in range(m+1)] for i in range(n+1)]

    # compute sum-of-pairs for each column of alignment beforehand to save time later
    sop = [0 for i in range(n)]
    for iKey in range(0, len(keys)-1):
        s = aln[keys[iKey]]
        for jKey in range(iKey+1, len(keys)):
            t = aln[keys[jKey]]
            for i in range(n):
                if (s[i] == '-' and t[i] != '-') or (s[i] != '-' and t[i] == '-'):
                    sop[i] += gap
                elif s[i] != '-' and t[i] != '-':
                    sop[i] += M[s[i]][t[i]]

    # perform DP
    for i in range(1,n+1):
        S[i][0] = S[i-1][0] + sop[i-1] + len(keys)*gap
        B[i][0] = 'i'
    for j in range(1,m+1):
        S[0][j] = S[0][j-1] + len(keys)*gap
        B[0][j] = 'j'
    for i in range(1,n+1):
        for j in range(1,m+1):
            options = {'i': S[i-1][j] + sop[i-1] + len(keys)*gap, 'j': S[i][j-1] + len(keys)*gap, 'ij': S[i-1][j-1] + sop[i-1]}
            for s in aln.values():
                if s[i-1] == '-':
                    options['ij'] += gap
                else:
                    options['ij'] += M[s[i-1]][seq[j-1]]
            S[i][j] = max(options.values())
            for o in options:
                if options[o] == S[i][j]:
                    B[i][j] = o
                    break

    # backtrack
    i,j = n,m
    out = {}
    for key in keys:
        out[key] = ''
    out[ID] = ''
    outScore = S[i][j]
    while i+j != 0:
        bt = B[i][j]
        if 'i' in bt:
            for key in keys:
                out[key] += aln[key][i-1]
            i -= 1
        else:
            for key in keys:
                out[key] += '-'
        if 'j' in bt:
            out[ID] += seq[j-1]
            j -= 1
        else:
            out[ID] += '-'
    for key in out:
        aln[key] = out[key][::-1]
    return outScore

# perform multiple sequence alignment without tree
def multiAlign(seqs, gap, treefile):
    if len(seqs) < 3:
        print("ERROR: Input contained less than 3 sequences")
        exit(-1)
    elif len(seqs) == 3:
        return tripleAlign(seqs, gap)

    # create distance matrix
    IDs = sorted(list(seqs.keys()))
    dm = {}
    for ID in IDs:
        dm[ID] = {ID:0}

    # generate distance matrix from tree
    if treefile != None:
        print("Tree file was specified. Computing distance matrix from tree...")
        try:
            import dendropy
        except:
            print("ERROR: Failed to import dendropy")
            print("You need dendropy installed to pass in a tree file")
            exit(-1)
        print("Reading tree...")
        tree = dendropy.Tree.get(file=treefile, schema='newick')
        print("Computing distance matrix...")
        pdm = tree.phylogenetic_distance_matrix()
        label2taxon = {}
        for leaf in tree.leaf_node_iter():
            label2taxon[leaf.taxon.label] = leaf.taxon
        for i in range(0, len(IDs)-1):
            for j in range(i+1, len(IDs)):
                iID = IDs[i]
                jID = IDs[j]
                d = pdm.patristic_distance(label2taxon[iID.replace('_',' ')], label2taxon[jID.replace('_',' ')])
                dm[iID][jID] = d
                dm[jID][iID] = d

    # generate distance matrix from pairwise alignments
    else:
        print("No tree file was specified. Computing distance matrix from pairwise alignments...")
        for i in range(0, len(IDs)-1):
            for j in range(i+1, len(IDs)):
                iID = IDs[i]
                jID = IDs[j]
                d = -1*pairAlign({iID:seqs[iID], jID:seqs[jID]}, gap)[0]
                dm[iID][jID] = d
                dm[jID][iID] = d

    # find 3 closest sequences for start
    print("Finding three closest sequences...")
    bestPair = None
    for i in range(0, len(IDs)-1):
        iID = IDs[i]
        for j in range(i+1, len(IDs)):
            jID = IDs[j]
            d = dm[iID][jID]
            if bestPair == None or d < bestPair[0]:
                bestPair = (d, iID, jID)
    iID,jID = bestPair[1],bestPair[2]
    bestThird = None
    for ID in IDs:
        if ID != iID and ID != jID:
            d = dm[iID][ID] + dm[jID][ID]
        else:
            d = float('inf')
        if bestThird == None or d < bestThird[0]:
            bestThird = (d, ID)
    kID = bestThird[1]

    # merge alignments
    print("Performing triple sequence alignment on three closest sequences...")
    out = tripleAlign({iID:seqs[iID], jID:seqs[jID], kID:seqs[kID]}, gap)[1]
    pq = Q.PriorityQueue()
    for ID in IDs:
        if ID not in out:
            pq.put((dm[iID][ID],ID))
            pq.put((dm[jID][ID],ID))
            pq.put((dm[kID][ID],ID))
    score = None
    print("Merging remaining sequences into alignment...")
    while len(out.keys()) < len(IDs):
        d,ID = pq.get()
        if ID in out:
            continue
        score = mergeAlign(out,ID,seqs[ID],gap)
        for newID in IDs:
            if newID not in out:
                pq.put((dm[newID][ID],newID))
    return score,out

# main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',   required=True,  type=argparse.FileType('r'),               help="Input FASTA File")
    parser.add_argument('-o', '--output',  required=True,  type=argparse.FileType('w'),               help="Output FASTA File")
    parser.add_argument('-p', '--indel',   required=False, type=float,                  default=0.1,  help="Indel Ratio")
    parser.add_argument('-t', '--tree',    required=False, type=argparse.FileType('r'), default=None, help="Guide Tree")
    args = parser.parse_args()
    start_time = time.time()
    print("User-specified indel rate: " + str(args.indel))
    gap = indelToGap(args.indel)
    print("Indel rate converted to gap rate: " + str(gap))
    print("Parsing input FASTA file...")
    seqs = parseFASTA(args.input)
    print("Computing performing sequence alignment...")
    aln = multiAlign(seqs, gap, args.tree)[1]
    print("Writing output alignment...")
    for ID in sorted(aln.keys()):
        args.output.write('>')
        args.output.write(ID)
        args.output.write('\n')
        args.output.write(aln[ID])
        args.output.write('\n')
    print("--- Completed in %s seconds ---" % (time.time() - start_time))
