import re 
import sys,os
import constant

# 拷贝config文件，并修改参数
# 拷贝模型文件
# 拷贝筛选特征

def readFileWriteTo(read_path, write_path):
    read_f = open(read_path,"rb")
    write_f = open(write_path, "wb")
    c = read_f.read()
    write_f.write(c)


def copyFeat(user_hash):
    # 拷贝筛选特征
    templete_dir = constant.TEMP_FEAT_DIR
    user_root_dir = constant.USER_PRED_ROOT_DIR
    for file in os.listdir(templete_dir):
        file_from_path = "%s/%s"%(templete_dir, file)
        file_to_path = "%s/%s/pipeline/%s"%(user_root_dir, user_hash, file)
        readFileWriteTo(file_from_path, file_to_path)

def copyModel(user_hash):
    # 拷贝模型
    model_from_path = constant.TEMP_MODEL_PATH
    user_root_dir = constant.USER_PRED_ROOT_DIR
    model_to_path = "%s/%s/result/ihalo_model.pickle"%(user_root_dir, user_hash)
    readFileWriteTo(model_from_path, model_to_path)

def copyConfig(user_hash):
    # 拷贝运行配置
    config_from_path = constant.TEMP_CONFIG_PATH
    user_root_dir = constant.USER_PRED_ROOT_DIR
    config_to_path = "%s/%s/result/pred_config.yml"%(user_root_dir, user_hash)
    readFileWriteTo(config_from_path, config_to_path)
    # 修改配置
    f = open(config_to_path,encoding="UTF8")
    c = f.read()
    new_c = re.sub("user_hash",user_hash,c)
    f.close()
    f  = open(config_to_path,"w",encoding="UTF8")
    f.write(new_c)

def copyCoreFeat(user_hash):
    # 拷贝核心特征
    from_path = constant.TEMP_CORE_FEAT_PATH
    user_root_dir = constant.USER_PRED_ROOT_DIR
    to_path = "%s/%s/result/core_select.csv"%(user_root_dir, user_hash)
    readFileWriteTo(from_path, to_path)

def main():
    user_hash = sys.argv[1]
    copyFeat(user_hash)
    copyModel(user_hash)
    copyConfig(user_hash)
    copyCoreFeat(user_hash)

if __name__=="__main__":
    main()