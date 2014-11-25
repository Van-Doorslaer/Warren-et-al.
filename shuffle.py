import random, re, itertools
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np

DNA=["a","c","t","g"]
motifs=[]
for r in range(2,4):
    iter= itertools.product(DNA,repeat = r)
    for i in iter:
        motifs.append("".join(i))




with open("FW.sliding.csv", "w") as out:
    print >>out, "virus, motif, window, under, over, neutral"
    for seq_record in SeqIO.parse("WG.fas","fasta"):
        print seq_record.id
        #run sliding window, save subsequences in "slides" list
        slides=[]
        window = 1000
        overlap = 100
        for x in range(0, len(seq_record.seq)+overlap, overlap):
            if x+window < len(seq_record.seq):
                slides.append( seq_record.seq[x:x+window] )
                x=x+overlap
            else:
                slides.append( seq_record.seq[x:]+seq_record.seq[:window-len(seq_record.seq[x:])] )
                x=x+overlap
                
        #randomize all subsequences 1000 times. Save the original + random in a list of lists (all_seq)
        all_seq=[]
        for record in slides:
            shuffled_recs=[record]
            nuc_list = list(record)
            for i in range(1000):
                random.shuffle(nuc_list)
                shuffled_recs.append( "".join(nuc_list) )
            all_seq.append(shuffled_recs)
        
        #count occurence of each motif in each original subsequence and in its random versions. 
        for m in motifs:
            n=0
            under=0
            over=0
            neutral=0
            
            for slice_set in all_seq:
                shuffled_count=[]
                not_shuffled_count=len(re.findall(m,str(slice_set[0]))) 
                for sr in slice_set[1:]:
                    shuffled_count.append(len(re.findall(m,str(sr))))
                a=np.array(shuffled_count)
    
                if not_shuffled_count <= np.percentile(a, 5):
                    under=under+1
                    n=n+1
                elif not_shuffled_count >= np.percentile(a, 95):
                    over=over+1
                    n=n+1
                else:
                    neutral=neutral+1
                    n=n+1
                
            results=(seq_record.id.split("|")[3].split("REF")[0], m,str(n), str(100*(float(under)/float(n))), str(100*(float(over)/float(n))), str(100*(float(neutral)/float(n))))
            print >>out, ",".join(results)
            #print ",".join(results)


