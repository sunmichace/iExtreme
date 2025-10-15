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
def load_json(file_path):
    # test "rhea_reaction_json/RHEA_10136.json"
    f = open(file_path, "rb")
    c = f.read().decode("UTF-8")
    return json.loads(c)

def get_vf(infile):
    f = open(infile)
    f.readline()
    table = []
    for line in f:
        arr = line.split("\t")
        table.append(arr[0])
    return table 

def to_vf_name(vf_map, vf_ids):
    table = []
    for vf_id in vf_ids:
        vf_name = vf_map[vf_id]
        table.append(vf_name)
    return table

def start_count(vf_names):
    table = {}
    for line in vf_names:
        strain, prot = line.split(" | ")
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
    map_table_path = "/home/LeiWei/bio_model/MetaVF_toolkit/MetaVF_toolkit-main/databases/vf_info_file.json"
    vf_map = load_json(map_table_path)
    vf_ids = get_vf(infile)
    vf_names = to_vf_name(vf_map, vf_ids)
    vf_count = start_count(vf_names)
    format_out(vf_count, outfile)

# conda run -n MetaVF_toolkit python metaVF.py -p $(pwd) -pjn draft_test -id example/E.coli_cds -outd test_draft_cds  -m draft -c 10 -ti 90 -tc 80

if __name__=="__main__":
    main()