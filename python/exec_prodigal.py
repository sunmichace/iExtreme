import os,sys 

# Test env in COMDEL2

# Example: python exec_prodigal.py GCA_001563815.1_ASM156381v1_genomic.fna prodigal_out

def exec(genome, outdir):
    command = "conda run -n fastani_prodigal prodigal -a %s/user_pred.pep -d %s/user_pred.cds -f gff -g 11 -o %s/user_pred.gff -p single -s %s/user_pred.stat -c -m -i %s 1>/dev/null 2>/dev/null"%(outdir, outdir, outdir, outdir, genome)
    os.system(command)

def main():
    user_genome = sys.argv[1]
    out_dir = sys.argv[2]
    exec(user_genome, out_dir)

if __name__ == "__main__":
    main()

# prokka --locustag frag --force --centre my --outdir prokka_out --prefix user_pred --kingdom Bacteria --cpus 1 GCA >/dev/null