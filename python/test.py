

def parseFastANI(text):
    text_str = text.split("fastANI res: ")[-1]
    name,ident,cov = text_str.split(",")
    return name,ident,cov

def parseVF(text):
    # VF res: 138
    num = int(text.split("VF res: ")[-1])
    return num

def parseARGs(text):
    # ARGs res: 2
    num = int(text.split("ARGs res: ")[-1])
    return num

def parsePGs(text):
    # PG res: 368
    num = int(text.split("PG res: ")[-1])
    return num

def parsePPRs(vf,pg,args):
    # PPRs
    num = (vf**2+pg**2+args**2)**0.5
    return num

print(parsePPRs(2,138,368))