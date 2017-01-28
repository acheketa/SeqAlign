#! /usr/bin/env python
'''
Triple sequence alignment using BLOSUM62

Niema Moshiri 2017
'''
from common import BLOSUM62 as M
from common import indelToGap
from common import parseFASTA
import argparse,sys,time

# align three sequences
def tripleAlign(seqs, g):
    # set up
    if len(seqs) != 3:
        sys.stderr.write("ERROR: Input did not contain exactly 3 sequences\n")
        exit(-1)
    rID,sID,tID = seqs.keys()
    r,s,t = seqs[rID],seqs[sID],seqs[tID]
    S = [[[0    for z in range(len(t)+1)] for y in range(len(s)+1)] for x in range(len(r)+1)]
    B = [[[None for z in range(len(t)+1)] for y in range(len(s)+1)] for x in range(len(r)+1)]

    # 1D base cases
    for x in range(1, len(r)+1):
        S[x][0][0] = 2*x*g
        B[x][0][0] = 'x'
    for y in range(1, len(s)+1):
        S[0][y][0] = 2*y*g
        B[0][y][0] = 'y'
    for z in range(1, len(t)+1):
        S[0][0][z] = 2*z*g
        B[0][0][z] = 'z'

    # 2D recursions
    for x in range(1, len(r)+1):
        for y in range(1, len(s)+1):
            options = {'xy': S[x-1][y-1][0] + M[r[x-1]][s[y-1]] + 2*g, 'x': S[x-1][y][0] + 2*g, 'y': S[x][y-1][0] + 2*g}
            S[x][y][0] = max(options.values())
            for o in options:
                if options[o] == S[x][y][0]:
                    B[x][y][0] = o
                    break
    for x in range(1, len(r)+1):
        for z in range(1, len(t)+1):
            options = {'xz': S[x-1][0][z-1] + M[r[x-1]][t[z-1]] + 2*g, 'x': S[x-1][0][z] + 2*g, 'z': S[x][0][z-1] + 2*g}
            S[x][0][z] = max(options.values())
            for o in options:
                if options[o] == S[x][0][z]:
                    B[x][0][z] = o
                    break
    for y in range(1, len(s)+1):
        for z in range(1, len(t)+1):
            options = {'yz': S[0][y-1][z-1] + M[s[y-1]][t[z-1]] + 2*g, 'y': S[0][y-1][z] + 2*g, 'Z': S[0][y][z-1] + 2*g}
            S[0][y][z] = max(options.values())
            for o in options:
                if options[o] == S[0][y][z]:
                    B[0][y][z] = o
                    break

    # 3D recursion
    for x in range(1, len(r)+1):
        for y in range(1, len(s)+1):
            for z in range(1, len(t)+1):
                options = {'xyz': S[x-1][y-1][z-1] + M[r[x-1]][s[y-1]] + M[r[x-1]][t[z-1]] + M[s[y-1]][t[z-1]],
                            'xy': S[x-1][y-1][z] + M[r[x-1]][s[y-1]] + 2*g, 'xz': S[x-1][y][z-1] + M[r[x-1]][t[z-1]] + 2*g, 'yz': S[x][y-1][z-1] + M[s[y-1]][t[z-1]] + 2*g,
                            'x': S[x-1][y][z] + 2*g, 'y': S[x][y-1][z] + 2*g, 'z': S[x][y][z-1]}
                S[x][y][z] = max(options.values())
                for o in options:
                    if options[o] == S[x][y][z]:
                        B[x][y][z] = o
                        break

    # backtrack (create strings in reverse order, and reverse them before returning)
    out = {rID:'', sID:'', tID:''}
    x,y,z = len(r),len(s),len(t)
    outScore = S[x][y][z]
    while sum((x,y,z)) != 0:
        bt = B[x][y][z]
        if 'x' in bt:
            out[rID] += r[x-1]
            x -= 1
        else:
            out[rID] += '-'
        if 'y' in bt:
            out[sID] += s[y-1]
            y -= 1
        else:
            out[sID] += '-'
        if 'z' in bt:
            out[tID] += t[z-1]
            z -= 1
        else:
            out[tID] += '-'
    for ID in out:
        out[ID] = out[ID][::-1]
    return outScore, out

# main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',  required=False, type=argparse.FileType('r'), default=sys.stdin,  help="Input FASTA File")
    parser.add_argument('-o', '--output', required=False, type=argparse.FileType('w'), default=sys.stdout, help="Output FASTA File")
    parser.add_argument('-p', '--indel',  required=False, type=float,                  default=0.1,        help="Indel Ratio")
    args = parser.parse_args()
    start_time = time.time()
    sys.stderr.write("User-specified indel rate: " + str(args.indel) + '\n')
    gap = indelToGap(args.indel)
    sys.stderr.write("Indel rate converted to gap rate: " + str(gap) + '\n')
    sys.stderr.write("Parsing input FASTA file...\n")
    seqs = parseFASTA(args.input)
    sys.stderr.write("Computing 3D dynamic programming matrix...\n")
    score,aln = tripleAlign(seqs, gap)
    sys.stderr.write("Writing output alignment...\n")
    for ID in sorted(aln.keys()):
        args.output.write('>')
        args.output.write(ID)
        args.output.write('\n')
        args.output.write(aln[ID])
        args.output.write('\n')
    sys.stderr.write("--- Completed in %s seconds ---" % (time.time() - start_time) + '\n')