import re
import sys

# 检验上传的文件是否合规
# 仅包含ATCG

def isGenome(path):
    try:
        f = open(path)
        f.readline()
    except:
        return False
    line = f.readline().replace("\"","")
    line = line.strip()
    line = line.split("\t")
    seq = line[1]
    if re.search("error",seq):
        return False
    return True

def main():
    input_path = sys.argv[1]
    if isGenome(input_path):
        print("genome")
    else:
        print("not genome")

if __name__=="__main__":
    main()
    # D:/阿伟的工作/7.酶查找数据库构建/1.数据收集/src/web/user_pred/a02308ba0bf006affbc24b09109b5d3c/raw/user_process.csv