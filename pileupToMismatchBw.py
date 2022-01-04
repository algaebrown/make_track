import sys
import pyBigWig
import pandas as pd


def count_mismatch(string):
    return sum([string.count(x) if type(x)==str else 0 for x in 'ATCGatcg'])

if __name__=='__main__':
    pileupf=sys.argv[1]
    chrsize=sys.argv[2]
    out=sys.argv[3]

    # read file in chunks
    chunks = pd.read_csv(pileupf, sep = '\t', header = None, names = ['chrom', 'pos', 'refbase', 'coverage', 'bases', 'a', 'b'], chunksize = 10**6)

    bw = pyBigWig.open(out, "w")

    # write header 
    chrsize_df = pd.read_csv(chrsize, sep = '\t', header = None)
    bw.addHeader([tuple(x) for x in list(chrsize_df.to_records(index = False))])

    for chunk in chunks:
        
        chunk['mismatch']=chunk['bases'].apply(count_mismatch)
        chunk = chunk.loc[chunk['mismatch']>0]
        bw.addEntries(chunk['chrom'].tolist(), chunk['pos'].tolist(), ends=(chunk['pos']+1).tolist(), values=chunk['mismatch'].astype(float).tolist())
    bw.close()
