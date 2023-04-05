import sys

if len(sys.argv) != 4:
    print(f'Usage: {sys.argv[0]} <input fasta file> <input AED file> <max AED score>', file=sys.stderr) 
    sys.exit(1)

scores = {}
score = float(sys.argv[3])

# read all of the AED scores calculated with get_names_AED.py
with open(sys.argv[2], 'r') as in_AED:
    for line in in_AED:
        line = line.strip()
        key, value = line.split()
        scores[key] = float(value)

write = False

with open(sys.argv[1], 'r') as in_fasta:
    for line in in_fasta:
        line = line.strip()

        # Header line
        if line[0] == '>':

            # get fasta id
            id = line[1:].split()[0]

            # if the id has an AED score better than the threshold, all subsequent lines should be written to stdout, otherwise not
            if id not in scores.keys() or scores[id] > score:
                write = False
            else:
                write=True
                
        # write the current line if it belongs to an entry with a good AED score
        if write:
            print(line)

    
