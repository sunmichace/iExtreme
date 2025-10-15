from sklearn.metrics import classification_report
from sklearn.model_selection import cross_validate, GridSearchCV
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV
from sklearn import svm
import sklearn
import numpy as np
import pandas as pd
import pickle
from time import sleep
import argparse,os

def printWaitMes(mes):
    print(mes)
    sleep(3)

class CodonSVM:
    def __init__(self, model_name, optimal_label,test_rate = 0.2):
        self.model_name = model_name
        self.model = svm.SVC()
        self.optimal_label = optimal_label
        self.optimal_model = None
        self.mean = None
        self.var = None 
        self.isstandard = False
        self.seed_pool = self.loadSeedPool()
        self.test_rate = test_rate
        self.model_save_path = 'optimal/%s.pickle'%(self.model_name)
        self.param = {
            'kernel':  ["rbf","poly","sigmoid"],#
            'C':  [1, 5, 10, 20, 40, 80,120],
            "gamma" : [0, 0.0001, 0.001, 0.1, 1, 10]
            #'epsilon': [0,0.01,0.1,1],
        }
        self.scoring = {
            "my_rule": "accuracy"
        }
        self.seed = 6294 # 随机取种子
        self.optimal_seed = None
        self.model_grid = GridSearchCV(self.model, self.param, cv=10, n_jobs=20,
                                        refit="my_rule", scoring=self.scoring)

    def saveOptimalModel(self,optimal_data):
        with open(self.model_save_path, 'wb') as f:  
            pickle.dump(optimal_data, f)
    
    def loadOptimalModel(self):
        with open(self.model_save_path, 'rb') as f:  
            optimal_model = pickle.load(f) 
            return optimal_model

    def standard(self,x):
        if not self.isstandard:
            self.mean = np.mean(x,axis=0)
            self.var = np.var(x,axis=0)
            self.isstandard = True
        return (x-self.mean)/self.var**0.5

    def loadSeedPool(self):
        f = open("seed_pool.txt","r")
        pool = f.readlines()
        pool = [int(seed.replace("\n","")) for seed in pool]
        return pool
    
    def paramOptimal(self,x,y):
        """
        param optimal
        """
        # split data
        std_x = self.standard(x)
        x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(std_x,y,test_size=self.test_rate, random_state = self.seed)
        # optimal start
        print("> Start to train and optimalize model......")
        res_grid = self.model_grid.fit(x_train, y_train)
        # save
        grid_record = res_grid.cv_results_
        best_index = res_grid.best_index_
        grid_record["my_best_index"] = best_index
        best_model = res_grid.best_estimator_
        optimal_model = {"model":best_model, "grid":grid_record,"mean":self.mean,"var":self.var}
        
        return optimal_model

    def optimalSeed(self,x,y,best_param):
        """
        random seed opeimal
        """
        std_x = self.standard(x)
        f1_score_list = []
        print("------------------------")
        print("index\tseed\tf1")
        print("------------------------")
        cur_i = 0
        for seed in self.seed_pool:
            cur_i += 1
            model = svm.SVC(
                C=best_param["C"],
                gamma=best_param["gamma"],
                kernel=best_param["kernel"]
            )
            x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(std_x,y,test_size=self.test_rate, random_state = seed)
            model.fit(x_train,y_train)
            seed_pred = model.predict(x_test)
            seed_report = classification_report(y_test, seed_pred,output_dict=True)
            f1_score = seed_report[self.optimal_label]["f1-score"]
            f1_score_list.append(f1_score)
            print("\r                            ",end="")
            print("\r%s\t%s\t%.3f"%(cur_i,seed,f1_score),end="")
        print()
        seed_df = pd.DataFrame({"seed":self.seed_pool,"f1":f1_score_list})
        seed_df = seed_df.sort_values(by='f1', ascending=False)
        length = seed_df.shape[0]
        median_i = length//2
        optimal_seed = seed_df.iloc[median_i,0]
        optimal_f1 = seed_df.iloc[median_i,1]
        print("> Optimal seed=%s, optiaml f1=%.3f"%(optimal_seed, optimal_f1))
        print("> Max f1=%.3f, min f1=%.3f, avearge f1=%.3f"%(seed_df["f1"].max(), seed_df["f1"].min(),seed_df["f1"].mean()))
        self.optimal_seed = optimal_seed
        return optimal_seed

    def finalModelTest(self,x,y,best_param,best_seed):
        """
        final test
        """
        std_x = self.standard(x)
        x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(std_x,y,test_size=self.test_rate, random_state = best_seed)
        model = svm.SVC(
                probability=True,
                C=best_param["C"],
                gamma=best_param["gamma"],
                kernel=best_param["kernel"]
        )
        model.fit(x_train,y_train)
        pred_train = model.predict(x_train)
        pred_test = model.predict(x_test)
        report_train = classification_report(y_train, pred_train,output_dict=False)
        report_test = classification_report(y_test, pred_test,output_dict=False)

        print("--"*10,"training data report","--"*10)
        print(report_train)
        sleep(3)
        # test
        print("--"*10,"test data report","--"*10)
        print(report_test)
        return model

    def train(self,x,y):
        """
        step1: param optimal
        step2: seed optimal
        step3: final model
        """
        #--------------------------------
        # step1: param optimal
        optimal_param = self.paramOptimal(x,y)
        best_param_index = optimal_param["grid"]["my_best_index"]
        best_param = optimal_param["grid"]["params"][best_param_index]
        print("> Model param optimization complete")
        sleep(3)
        print("> Best param is",best_param)
        sleep(3)

        #-------------------------
        # step2: seed optimal 
        print("> Start to optimal seed[%d] for label %s......"%(len(self.seed_pool), self.optimal_label))
        sleep(3)
        best_seed = self.optimalSeed(x,y,best_param)

        #--------------------------------
        # step3 
        print("> Final model test......")
        sleep(3)
        final_model = self.finalModelTest(x,y,best_param,best_seed)
        # save
        optimal_data = {"model":final_model, "optimal_history":optimal_param,"mean":self.mean,"var":self.var}
        self.saveOptimalModel(optimal_data)
        print("> Final model save to",self.model_save_path)

    def predict(self,x,y=None,type="invalid"):
        """
        predict
        """
        # load
        optimal_model = self.loadOptimalModel()

        # standard
        std_x = (x-optimal_model["mean"])/optimal_model["var"]**0.5
        pred_y = optimal_model["model"].predict(std_x)
        pred_y_proba = optimal_model["model"].predict_proba(std_x)
        if type=="valid":
            report = classification_report(pred_y, y,output_dict=False)
            print("--"*10,"predict data report","--"*10)
            print(report)
        return pred_y, pred_y_proba

def getArg():
    parser = argparse.ArgumentParser(description='codon+SVM model to predict extreme microorganism')
    parser.add_argument('--model', type=str, required=True,help='model name')  # , choices=["exheat","exsalt","expH","all"]
    parser.add_argument('--model_type', required=True,choices=["train","predict","predict_all"], help='train or predict')
    parser.add_argument('--data', required=True,type=str,help='data location must be csv format')
    parser.add_argument('--optimal_label', required=True, default="1",choices=["0","1","2","all"],help='optimal label')
    parser.add_argument('--out', default="empty",type=str,help='out predict result')
    parser.add_argument('--test_rate', default=0.2,type=float,help='test data rate')
    parser.add_argument('--proba', default=True,action='store_true',help='out probability')
    args = parser.parse_args()
    return args

def train(args):
    model_name = args.model
    optiaml_label = args.optimal_label
    test_rate = args.test_rate
    # --------------------------
    codon_svm = CodonSVM(model_name, optiaml_label, test_rate)
    df = pd.read_csv(args.data,index_col=0,sep="\t")
    X = df.iloc[:,0:-1].values
    Y = df.iloc[:,-1].values
    codon_svm.train(X,Y)

def predict(args):
    model_name = args.model
    optiaml_label = args.optimal_label
    test_rate = args.test_rate
    proba = args.proba
    out = args.out
    assert out!="empty","predict mode must point --out"
    # --------------------------
    codon_svm = CodonSVM(model_name, optiaml_label, test_rate)
    df = pd.read_csv(args.data,index_col=0,sep="\t")
    X = df.iloc[:,0:-1].values
    Y = df.iloc[:,-1].values
    print("> Predict start......")
    label, label_proba = codon_svm.predict(X,type="invalid")
    ex_index = int(optiaml_label)
    tmp = {"name":df.index} # ,"ex_proba":ex_proba
    for index in range(len(label_proba[0])): # [ex_index]
        key_name = "proba_label_%d"%(index)
        tmp[key_name] = label_proba[:,index]
    pred_df = pd.DataFrame(tmp)
    pred_df.to_csv(out,index=None)
    print("> Predict complete")
    print("> Predict data save to ",out,end="\n\n")

def predict_all(args):
    model_name = args.model
    optiaml_label = args.optimal_label
    test_rate = args.test_rate
    proba = args.proba
    out = args.out
    assert out!="empty","predict mode must point --out"
    # --------------------------
    codon_svm_exheat = CodonSVM("exheat", "2", test_rate)
    codon_svm_exsalt = CodonSVM("exsalt", "1", test_rate)
    codon_svm_expH = CodonSVM("expH", "1", test_rate)
    df = pd.read_csv(args.data,index_col=0,sep="\t")
    X = df.iloc[:,0:-1].values
    Y = df.iloc[:,-1].values
    print("> Predict start......")
    label_exheat, label_proba_exheat = codon_svm_exheat.predict(X,type="invalid")
    label_exsalt, label_proba_exsalt = codon_svm_exsalt.predict(X,type="invalid")
    label_expH, label_proba_expH = codon_svm_expH.predict(X,type="invalid")
    # 
    ex_proba_heat = label_proba_exheat[:,2]
    ex_proba_salt = label_proba_exsalt[:,1]
    ex_proba_pH = label_proba_expH[:,1]
    pred_df = pd.DataFrame({
        "name":df.index,
        "exheat_proba":ex_proba_heat,
        "exsalt_proba":ex_proba_salt,
        "expH_proba":ex_proba_pH,
    })
    pred_df.to_csv(out,index=None)
    print("> Predict complete")
    print("> Predict data save to ",out,end="\n\n")

def main():
    args = getArg()

    model_files_path = {
        "expH":"data/ccdon_freq_selected_pH_less5.tsv",
        "exsalt":"data/ccdon_freq_selected_salinity.tsv",
        "exheat":"data/ccdon_freq_selected_temp_large70.tsv",
    }
    assert os.path.exists(args.data), "%s not exists"%(args.data)
    if args.model_type=="train":
        train(args) 
    elif args.model_type=="predict":
        predict(args) 
    elif args.model_type=="predict_all":
        predict_all(args) 

if __name__=="__main__":
    main()

"""
python codon_svm.py \
    --model iMicrobes \
    --model_type predict \
    --data example/test.tsv \
    --optimal_label 2 \
    --out example/predict.csv
"""