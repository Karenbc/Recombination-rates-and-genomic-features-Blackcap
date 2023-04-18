#!/usr/bin/env python3
import sys

# seperator of the output, change to 

if len(sys.argv) != 3:
    print(f'Usage: {sys.argv[0]} <recombination file> <window size>')
    sys.exit(1)

try:
    window_size = int(sys.argv[2])
except ValueError:
    print(f'Usage: {sys.argv[0]} <recombination file> <window size>')
    sys.exit(1)

curr_window = []
curr_chr = ''

window_start = 0

with open(sys.argv[1], 'r') as in_file:
    # go through the file
    for line in in_file:
        line = line.strip()
        data = line.split()

        # get data from line
        chr, begin, end, rate = data[0], int(data[1]), int(data[2]), float(data[3])

        if curr_chr != chr:
            if curr_chr != '':
                mean = sum([r * (e - b) for e, b, r in curr_window]) / sum([(e - b) for e, b, r in curr_window])
                print('\t'.join([curr_chr, str(curr_window[0][0]), str(curr_window[-1][1]), str(mean)]))
            curr_chr = chr
            curr_window = [(0, end, rate)]
            window_start = 0
        
        while end > window_start + window_size:
            curr_window.append((begin, window_start + window_size, rate))
            begin = window_start + window_size

            mean = sum([r * (e - b) for e, b, r in curr_window]) / sum([(e - b) for e, b, r in curr_window])
            print('\t'.join([curr_chr, str(curr_window[0][0]), str(curr_window[-1][1]), str(mean)]))

            window_start += window_size
            curr_window = []
            
        curr_window.append((begin, end, rate))


mean = sum([r * (e - b) for e, b, r in curr_window]) / sum([(e - b) for e, b, r in curr_window])
print('\t'.join([curr_chr, str(curr_window[0][0]), str(curr_window[-1][1]), str(mean)]))

        
