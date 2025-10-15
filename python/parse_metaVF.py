import sys,os,re,json

# 解析metaVF结果
"""
{results:[
    {
        name:E.coli,
        children:[
            {
                name:protein,
                value:12
            }
        ]
    },
    {},
]}
"""
def parse(infile):
    f = open(infile)
    f.readline()
    table = {}
    for line in f:
        arr = line.split("\t")
        strain, prot = arr[0].split(" | ")
        prot = re.sub(r"\(.+?\)", "", prot)
        prot = prot.strip(" ")
        if strain not in table:
            table[strain] = {prot:0}
        if prot not in table[strain]:
            table[strain][prot] = 0
        table[strain][prot] += 1
    return table 

def format_out(result, outfile):
    table = {"results":[]}
    for i in result:
        tmp = {"name":i, "children":[]}
        for j in result[i]:
            tmp["children"].append({"name":j, "value":result[i][j]})
        table["results"].append(tmp)
    if len(table["results"])==0:
        table["results"].append(
            {"name":"No VF","children":[
                {"name":"No VF Found","value":1}
            ]}
        )
    # save
    out_f = open(outfile, "w") 
    out_f.write(json.dumps(table))

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    result = parse(infile)
    format_out(result, outfile)

# conda run -n MetaVF_toolkit python metaVF.py -p $(pwd) -pjn draft_test -id example/E.coli_cds -outd test_draft_cds  -m draft -c 10 -ti 90 -tc 80

if __name__=="__main__":
    main()