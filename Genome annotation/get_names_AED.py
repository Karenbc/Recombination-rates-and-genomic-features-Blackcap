from os import sep
import sys

if len(sys.argv) != 2:
    print(f'Usage: {sys.argv[0]} <gff file>', file=sys.stderr)
    sys.exit(1)


with open(sys.argv[1], 'r') as in_gff:
    for line in in_gff:
        # ignore comment lines
        if line[0] == '#': 
            continue
        
        # gff is tab delimited
        line = line.strip()
        line = line.split('\t')

        # only the mRNA fields have AED scores
        if line[1] != 'maker' or line[2] != 'mRNA':
            continue

        # AED score and name is in the last column which contains semicolon delimited key-value pairs
        data = [d.split('=') for d in line[-1].split(';') if d != '']
        data = dict(data)

        # print the ID and the AED score
        print(data['ID'], data['_AED'], sep='\t')
