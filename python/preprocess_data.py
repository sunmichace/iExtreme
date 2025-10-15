import re
import sys

# process data
# format as follow
"""
col = {refseq, 4096 feats, label}
"""

def read_codon_table():
    f = open("feat.txt")
    feat = {feat.strip():0 for feat in f}
    return feat

def count_freq_codon(seqs):
    codon_table = read_codon_table()
    totol_codon = 0
    for seq in seqs:
        seq = re.sub(r"\s","", seq)
        for i in range(0,len(seq),3):
            codon_2 = seq[i:i+6]
            if len(codon_2)!=6:
                continue
            if codon_2 not in codon_table:
                continue
            totol_codon += 1
            codon_table[codon_2] += 1
    for k in codon_table:
        codon_table[k]/=totol_codon
    return codon_table

def getOneSeqs(path):
    f = open(path)
    try:
        c = f.read()
    except:
        c = ">error\nerror\n"
    m = re.findall(">.+\n([^>]+)", c)
    # seq = "".join(m)
    # seq = re.sub("\s","",seq)
    seqs = m
    return seqs

def outFile(name, ccodon_freq_table, out_path):
    f = open(out_path, "w")
    feat = list(ccodon_freq_table.keys())
    freq = list(ccodon_freq_table.values())
    freq = list(map(str, freq))
    head_name = "\t".join(feat)+"\t"+"label\n"
    freq_line = "\t".join(freq)+"\t"+"-1"
    f.write("refseq\t"+head_name)
    f.write("%s\t%s\n"%(name,freq_line))

def main():
    input_path = sys.argv[1]
    outfile_path = sys.argv[2]
    name = "user"
    seqs = getOneSeqs(input_path)
    if seqs[0] == "error\n":
        ccodon_freq_table = {"error":"error"}
    else:
        ccodon_freq_table = count_freq_codon(seqs)
    outFile(name, ccodon_freq_table, outfile_path)

if __name__=="__main__":
    main()