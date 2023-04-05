import sys

if len(sys.argv) != 4:
    print(f'Usage {sys.argv[0]} <in gff file> <in AED file> <max AED score>', file=sys.stderr)
    sys.exit(1)

score = float(sys.argv[3])

scores = {}

# read all of the AED scores calculated with get_names_AED.py
with open(sys.argv[2], 'r') as in_AED:
    for line in in_AED:
        line = line.strip()
        key, value = line.split()
        scores[key] = float(value)

exon_id = 1

print('##gff-version 3')

with open(sys.argv[1], 'r') as in_gff:
    for line in in_gff:
        line = line.strip()
        # ignore comment lines
        if line[0] == '#':
            continue

        # gff is tab delimited
        l = line.split('\t')

        # print scaffold lines
        if l[1] == '.': 
            print(line)
            continue
    
        # ignore non maker lines
        if l[1] != 'maker':
            continue

        # AED score and name is in the last column which contains semicolon delimited key-value pairs
        data = dict([d.split('=') for d in l[-1].split(';') if d != ''])

        # if the line is a gene, at least of the corresponding mRNAs must have a AED value lower than the threshold
        if l[2] == 'gene':
            i = 1
            exon_id = 1
            while True:
                # The mRNA ID is the Gene ID + mRNA-{number}
                key = f'{data["ID"]}-mRNA-{i}'
                # No more mRNAs corresponding to this gene
                if key not in scores.keys():
                    break
                # A mRNA was found that has a lower score than the threshold
                if scores[key] <= score:
                    print(line)
                    break
                i += 1

            continue

        # for mRNAs only the ID is important
        if l[2] == 'mRNA':
            if data['ID'] in scores and scores[data['ID']] <= score:
                print(line)
            continue

        # Features other than mRNA or gene (exon, CDS, ...) have a list of parents, one of which needs to have a correct score
        keys = data['Parent'].split(',')

        # Get all parents that have an acceptable score
        keys = [key for key in keys if key in scores.keys() and scores[key] <= score]

        # Replace parents with acceptable ones
        data['Parent'] = ','.join(keys)

        # rebuild line if feature still has parents 
        if len(keys) > 0:
            if l[2] == 'exon':
                data['ID'] = keys[0] + ':' + str(exon_id)

                l[-1] = ';'.join([f'{k}={v}' for k, v in data.items()]) + ';'
                line = '\t'.join(l)
                exon_id += 1

            print(line)