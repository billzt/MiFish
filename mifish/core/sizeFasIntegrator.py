from Bio.SeqIO.FastaIO import SimpleFastaParser

def run(zotusCountFile:str, zotusFaFile:str, resultFile:str):
    zotuCount = {}
    with open(zotusCountFile) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            (zotu, count) = line.strip().split()
            zotuCount[zotu] = count

    with open(zotusFaFile) as handle:
        with open(resultFile, mode='w') as out_handle:
            for (zotu, seq) in SimpleFastaParser(handle):
                if zotu in zotuCount:
                    size = zotuCount[zotu]
                    if size != '0':
                        print(f'>{zotu};size={zotuCount[zotu]};\n{seq}', file=out_handle)
