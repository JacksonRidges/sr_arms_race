#!/usr/bin/env python3

import sys

def main():
	print("position", "x", "n", "folded", sep="\t")
	for i, line in enumerate(sys.stdin):
		if i == 0:
			continue
		sl = line.rstrip('\n').split('\t')
		print("\t".join([sl[1], str(int(round(float(sl[5]) * float(sl[7])))), sl[7], "1"]))

if __name__ == "__main__":
	main()
