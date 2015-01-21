#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys
import argparse

class IndexedFasta(object):
	def __init__(self):
		self.count = 0
		self.dict = dict()

	def open(self, fileName):
		with open(fileName, 'r') as rep_set:
			if not rep_set:
				print("Could not open: ", rep_set, ".", file=sys.stderr)
				return False
			nextLine = rep_set.readline()
			if not nextLine:
				print("FASTA is empty.", file=sys.stderr)
				return False
			currentSequence = nextLine
			while True:
				nextLine = rep_set.readline()
				if nextLine.startswith('>') :
					#Start of the next sequence so process current one and start over
					self.processSequence(currentSequence)
					currentSequence = nextLine
					continue
				if not nextLine:
					#EOF so process the last sequence and break
					self.processSequence(currentSequence)
					break
				#Otherwise just concatinate to keep building the sequence
				currentSequence += nextLine.rstrip('\n')
			print("Indexed", self.count, "sequences...", file=sys.stderr)
		return True

	def processSequence(self, sequence):
		parts = sequence.split()
		if len(parts) < 2:
			print("Malformed FASTA entry:", sequence, file=sys.stderr)
			return
		otu = parts[0].translate(None, '>')
		rep_seq = parts[-1]
		self.dict[otu] = rep_seq
		self.count += 1

	def __getitem__(self, key):
		return self.dict[key]

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("rep_set", help="The representative sequence set in FASTA corresponding to the supplied OTU table.", type=str)
	parser.add_argument("otu_table", help="The OTU table to which to add the sequences. Must be in tab delimited format, not BIOM format.", type=str)
	args = parser.parse_args()

	#First we suck the FASTA file into memory
	print("Opening FASTA...", file=sys.stderr)
	indexedFasta = IndexedFasta()
	if not indexedFasta.open(args.rep_set):
		sys.exit(1)

	#Next we go through the OTU table and pair representative sequences to the OTUs
	print("Opening OTU table...", file=sys.stderr)
	otu_table = open(args.otu_table, 'r')
	if not otu_table:
		print("Could not open", args.rep_set, ".", file=sys.stderr)
		sys.exit(1)

	#Copy the first line exactly
	line = otu_table.readline()
	if not line:
		print("OTU table ended prematurely.", file=sys.stderr)
		sys.exit(1)
	print(line, end='', file=sys.stdout)

	#Add the new column description
	line = otu_table.readline()
	if not line:
		print("OTU table ended prematurely.", file=sys.stderr)
		sys.exit(1)
	print(line.rstrip('\n'), "\tReprSequence", sep='', file=sys.stdout)
	
	#Pair all the otu's to their representative sequences
	count = 0
	missed = 0 
	for line in otu_table:
		parts = line.split('\t')
		if len(parts) == 0:
			print("Malformed Line!", file=sys.stderr)
			continue
		#parts[0] is the otu
		try:
			sequence = indexedFasta[parts[0]]
			print(line.rstrip('\n'), '\t', sequence, sep='', file=sys.stdout)
			count += 1
			if count % 10000 == 0:
				print("Married", count, "sequences...", file=sys.stderr)
		except KeyError:
			print("OTU representative sequence not found in FASTA", parts[0], file=sys.stderr)
			print(line, end='', file=sys.stdout)
			missed += 1
	print("Married", count, "sequences...", file=sys.stderr)
	if missed:
		print("Unable to match", missed, "OTUs...", file=sys.stderr)


if __name__ == "__main__":
	main()
