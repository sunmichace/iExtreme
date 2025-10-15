import re,sys
import constant

# 解析预测的结果

"""预测文件内容
gene_id,label,proba
user,0.0,0.9905106583860737
"""

def parseStart(user_hash):
    # 返回是嗜盐菌的概率
    # USER_PRED_RESULT_PATH = "%s/../model/ihalo_user_pred_log/result/predict_labels.csv"%(os.getcwd())
    user_root_dir = constant.USER_PRED_ROOT_DIR
    pred_path = "%s/%s/result/user_predict.csv"%(user_root_dir, user_hash)
    f = open(pred_path)
    f.readline()
    line = f.readline().strip()
    arr = line.split(",")
    score_0 = float(arr[1])
    score_1 = float(arr[2])
    score_2 = float(arr[3])
    proba = -1*score_0+score_2
    return proba

def main():
    user_hash = sys.argv[1]
    proba = parseStart(user_hash)
    print(proba)

if __name__=="__main__":
    main()

