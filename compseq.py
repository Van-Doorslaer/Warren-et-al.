import itertools, os
from Bio import SeqIO



wordsize=[2,3]
for w in wordsize:
    dinuc=[]
    for i in itertools.product("ACGT", repeat=w):
        dinuc.append("".join(i))
    topline=[""]+dinuc

    #print dinuc
    print "calculating the frequency of sequences with wordsize "+str(w)
    with open(str(w)+"_subset_nucleotide_count.csv","w") as out:
        print >>out, ",".join(topline)
        for record in SeqIO.parse("subset.fas","fasta"):
            with open("temp.fasta.fas","w") as temp:
                ID=record.id.split("|")[3].split(".")[0]
                print >>temp, ">"+ID
                print >>temp, record.seq
            os.system( "compseq -sequence temp.fasta.fas -word "+str(w)+" -outfile temp.comp.txt -calcfreq Y" )
        
            with open("temp.comp.txt","rU") as f:
                results=[]
                results.append(ID)
                for line in f:
                    line=line.strip().split("\t")
                    if line[0] in dinuc:
                        results.append(line[5])
                        #print line
                print >>out, ",".join(results)
