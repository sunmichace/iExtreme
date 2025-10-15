import os,sys 

# Test env in muscle

# exec fastANI
# Example:python exec_fastani.py q_genome.fna r_genomes.txt ident cov_rate
#   q_genome.fna: query genome path
#   r_genomes.txt：file containing ref genomes path segement in \n
#   ident: ident threshold [0-100%]
#   cov_rate: cov_rate threshold [0-100%]

# out example: 
#   >Empty: fastANI res: NULL,NULL,NULL
#   >Normal: fastANI res: strain_name,75.5062,0.17366799921357888

def exec(q_genome, r_genomes, res):
    command = "fastANI -q %s --rl %s -o %s 1>/dev/null 2>/dev/null"%(q_genome, r_genomes, res)
    os.system(command)

def parse(res, ident_r, cov_rate_r):
    # get max ident, just first line
    f = open(res)
    first_line = f.readline().strip()
    q_g,r_g,ident,cov_num,all_num = first_line.split("\t")
    cov_rate = 100*float(cov_num)/float(all_num)
    r_gname = r_g.split("/")[-1]
    print_line = "fastANI res: NULL,NULL,NULL"
    if float(ident)>=ident_r and cov_rate>=cov_rate_r:
        print_line = "fastANI res: %s,%s,%s"%(r_gname, ident, cov_rate)
    return print_line

def main():
    q_genome = sys.argv[1]
    r_genomes = sys.argv[2]
    res = sys.argv[3]
    ident_threshold = float(sys.argv[4])
    cov_rate_threshold = float(sys.argv[5])
    exec(q_genome, r_genomes, res)
    final_res = parse(res, ident_threshold, cov_rate_threshold)
    print(final_res)

if __name__ == "__main__":
    main()