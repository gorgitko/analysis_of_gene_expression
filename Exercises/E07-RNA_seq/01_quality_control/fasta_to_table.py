#!/usr/bin/env python3

import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', '--delimiter', type=str, default='\t', help='CSV delimiter')
parser.add_argument('fasta_file', type=str, help='FASTA file')
args = parser.parse_args()

with open(args.fasta_file, mode='r') as f:
    rows = [x.strip() for x in f.readlines() if x.strip()]

for l1, l2 in zip(rows[::2], rows[1::2]):
    print('{}{}{}'.format(l1[1::], args.delimiter, l2))

print()
