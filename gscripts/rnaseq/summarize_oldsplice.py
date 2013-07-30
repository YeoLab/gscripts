import cPickle as pickle
import sys
exonsDetected = {}

minReadsTotal = 30
minReadsInEach =0

exonTypes = ["SE", "MXE"]
from collections import Counter, defaultdict
exonsDetected = defaultdict(Counter)

if __name__ == "__main__":
	f = pickle.load(open(sys.argv[1], 'rb'))
	for gene in f:
		thisExonTypes = set(gene.keys()).intersection(exonTypes)
		for exonType in thisExonTypes:
			for exon in gene[exonType].keys():
				thisLine = gene['descriptor'] + "\t" + exonType + "\t" + exon 
				thisReads = 0
				nIso = 0
				for isoform in sorted(gene[exonType][exon].keys()):
					if gene[exonType][exon][isoform] >= minReadsInEach:
						nIso += 1
					thisReads += gene[exonType][exon][isoform]
					thisLine += "\t" + str(gene[exonType][exon][isoform])
				print thisLine
				if thisReads >= minReadsTotal:
					exonsDetected[exonType]['detected'] += 1
				else:
					exonsDetected[exonType]['notDetected'] += 1
for ee in exonTypes:
	sys.stderr.write("\t".join(map(str, [ee, exonsDetected[ee]['detected'], exonsDetected[ee]['notDetected']])) + "\n" )	


			
