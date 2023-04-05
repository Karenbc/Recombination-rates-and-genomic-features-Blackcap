'''
This script relies on the QI values that maker prints in its out fasta files if the pred_stats option is set to 1.
The QI values have the following format:

QI:0|0.77|0.68|1|0.77|0.78|19|462|824
The fields stand for the following:
1. Length of the 5' UTR
2. Fraction of splice sites confirmed by an EST/mRNA-seq alignment
3. Fraction of exons that match an EST/mRNA-seq alignment
4. Fraction of exons that overlap EST/mRNA-seq or protein alignments
5. Fraction of splice sites confirmed by ab initio gene prediction
6. Fraction of exons that overlap an ab initio gene prediction
7. Number of exons in the mRNA
8. Length of the 3' UTR
9. Length of the protein sequence produced by the mRNA

The script only uses fields 1 and 8.
'''

import sys

if len(sys.argv) != 2:
    print(f'Usage: {sys.argv[0]} <maker generated fasta>')
    sys.exit(1)

# length of each line in the resulting fasta file
LINE_LENGTH=80


def write_seq(seq, header):
    'Writes a sequence in fasta format to stdout'

    # don't print an empty sequence
    if len(seq) == 0:
        return

    # print the header
    print(header)
    # print the sequence in lines
    for i in range(len(seq) // LINE_LENGTH + 1):
        print(seq[i * LINE_LENGTH:(i+1) * LINE_LENGTH])

seq = ''
begin = 0
end = 0
header = ''

# read the in fasta file
with open(sys.argv[1], 'r') as in_fasta:
    for line in in_fasta:

        # remove trailing whitespace and line end
        line = line.strip()

        # if the line is a header
        if line[0] == '>':
            # print the last sequence read
            seq = seq[begin:end]
            write_seq(seq, header)

            # reset the sequence
            seq = ''
            line = line.split(' ')
            # QI values are at the end of the fasta name
            data = line[-1][3:].split('|')
            # begin and end are at positions 0 and 7
            begin = int(data[0])
            # field 7 is the length of the 3' UTR -> use negative value for index from the end
            end = -int(data[7])
            # set length of UTRs to 0, as they are not in the output
            data[0] = '0'
            data[-2] = '0'
            # reassemble QI values
            line[-1] = 'QI:' + '|'.join(data)
            # remove offset value (there is no offset anymore)
            for i in range(len(line)):
                if line[i].split(':')[0] == 'offset':
                    break
            else:
                i = len(line)
            # reassemble header
            line = line[:i] + line[(i+1):]
            header = ' '.join(line)

        else:
            # no header -> sequence
            seq += line
