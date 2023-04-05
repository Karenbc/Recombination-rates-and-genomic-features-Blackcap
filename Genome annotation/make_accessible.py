import sys


if len(sys.argv) != 2:
    print(f'Usage {sys.argv[0]} <in_gff>', file=sys.stderr)
    sys.exit(1)

# print header
print('scaffold,start,end,strand,ID,function,interproID,PfamID,Ontology')

with open(sys.argv[1], 'r') as in_gff:
    for line in in_gff:
        line = line.strip()
        # don't print comments
        if line[0] == '#':
            continue
    
        # gff is tab delimited
        line = line.split('\t')

        # only print mRNAs
        if line[2] != 'mRNA':
            continue

        # annotations are in the last column which contains semicolon delimited key-value pairs
        data = dict([d.split('=') for d in line[-1].split(';') if d != ''])

        # check if homolog is known 
        f = data['Note']
        if f == 'Protein of unknown function':
            f = 'unknown function'
        else:
            # remove the 'Similar to ' from the homologs name
            f = f[len('Similar to '):]
        
        # There should be no commas in the gene names
        if ',' in f:
            print(f"I've found a comma in a protein name. This will disrupt the csv file. The protein is: {data['ID']}, please fix.", file=sys.stderr)

        # interpro and pfam domains can be multiple values
        ipro, pfam = [], []

        # interpro and pfam domains are kept in a combined field called "Dbxref" and are comma seperated in there
        if 'Dbxref' in data.keys():
            for key, value in [d.split(':') for d in data['Dbxref'].split(',')]:
                if key == 'InterPro' and value != '-':
                    ipro.append(value)
                if key == 'Pfam' and value != '-':
                    pfam.append(value)

        # Gene ontology can also consist of multiple terms
        ontology = []

        # The ontology terms are in a seperate field called "Ontology_term"
        if 'Ontology_term' in data.keys():
            for value in [d.split(':')[1] for d in data['Ontology_term'].split(',') if d != 'GO:-']:
                ontology.append(value)
        
        # print all the found values comma seperated and the fields with multiple results semi-colon seperated
        out = [line[0], line[3], line[4], line[6], data['ID'], f, ';'.join(ipro), ';'.join(pfam), ';'.join(ontology)]

        print(','.join(out))