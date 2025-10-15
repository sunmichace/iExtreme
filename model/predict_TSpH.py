import joblib 
import os,sys,re 
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from collections import Counter
import numpy as np
import pandas as pd 

def load_means_stds(predictor):
    """
    load means and stds from predictor file
    """
    means=dict()
    stds=dict()
    features=list()
    for line in open(predictor.replace('pkl','f'),'r'):
        if line.startswith('#'):continue
        cont=line.split()
        means[cont[0]]=float(cont[1])
        stds[cont[0]]=float(cont[2])
        features.append(cont[0])
    return means,stds,features

def load_model(path):
    """
    load model from path
    """
    model=joblib.load(path)
    means,stds,features = load_means_stds(path)
    return model,means,stds,features

def do_count(seq):
    """
    count the number of pep
    """
    dimers = Counter()
    for i in range(len(seq)-1): dimers[seq[i:i+2]] += 1.0
    return dimers

def count_dimer(fasta_file,p):
    seqs = [str(rec.seq).upper() for rec in SeqIO.parse(fasta_file,'fasta')]

    if int(p) == 0:num_cpus = cpu_count()
    else: num_cpus = int(p)
    results = Pool(num_cpus).map(do_count, seqs)
    dimers = sum(results, Counter())
    return dict(dimers)

def get_dimer_frequency(fasta_file,p):
    dimers = count_dimer(fasta_file,p)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
    dimers_fq = dict()
    for a1 in amino_acids:
        for a2 in amino_acids:
            dimers_fq[a1+a2] = dimers.get(a1+a2,0.0)
    number_of_aa_in_fasta = sum(dimers_fq.values())
    for key,value in dimers_fq.items(): dimers_fq[key] = value/number_of_aa_in_fasta
    return dimers_fq

def predict(fasta_file,model,means,stds,features,p=2):
    """
    user model to predict
    """
    dimers_fq = get_dimer_frequency(fasta_file,p)
    Xs = list()
    for fea in features:
        Xs.append((dimers_fq[fea]-means[fea])/stds[fea])
    Xs = np.array(Xs).reshape([1,len(Xs)])
    pred_ogt = model.predict(Xs)[0]
    return np.around(pred_ogt,decimals=2)

def get_diff_x_y(y, diff_path):
    """
    get std in confidence interval
    """
    df = pd.read_csv(diff_path)
    diff_x = df["diff_x"].values + y 
    diff_x[diff_x<0]=0
    diff_y = df["diff_p"].values
    return np.around(diff_x,decimals=2), np.around(diff_y,decimals=2)


def start(prot_fasta, out_file, models, diff_path):
    # load model
    temp_model,means_T,stds_T,features = load_model(models["Temp"])
    salinity_model,means_S,stds_S,features = load_model(models["Salinity"])
    pH_model,means_pH,stds_pH,features = load_model(models["pH"])
    # predict
    temp_pred = predict(prot_fasta,temp_model,means_T,stds_T,features)
    salinity_pred = predict(prot_fasta,salinity_model,means_S,stds_S,features)
    pH_pred = predict(prot_fasta,pH_model,means_pH,stds_pH,features)
    # add coef
    temp_diff_x, temp_diff_y = get_diff_x_y(temp_pred, diff_path["Temp"])
    salinity_diff_x, salinity_diff_y = get_diff_x_y(salinity_pred, diff_path["Salinity"])
    pH_diff_x, pH_diff_y = get_diff_x_y(pH_pred, diff_path["pH"])
    # print
    print(">Temp_x=",",".join(map(str,temp_diff_x)))
    print(">Temp_y=",",".join(map(str,temp_diff_y)))
    print(">Salinity_x=",",".join(map(str,salinity_diff_x)))
    print(">Salinity_y=",",".join(map(str,salinity_diff_y)))
    print(">pH_x=",",".join(map(str,pH_diff_x)))
    print(">pH_y=",",".join(map(str,pH_diff_y)))

def main():
    """
    main func
    """
    model_path = {
        "Salinity":"./optimal/train_salt.pkl",
        "Temp":"./optimal/train_temp.pkl",
        "pH":"./optimal/train_pH.pkl",
    }
    model_diff_path = {
        "Salinity":"./data/my_pred_train_df_salt.csv.diff",
        "Temp":"./data/my_pred_train_df_temp.csv.diff",
        "pH":"./data/my_pred_train_df_pH.csv.diff",
    }
    prot_fasta = sys.argv[1]
    out_file = sys.argv[2] 
    start(prot_fasta, out_file, model_path, model_diff_path)


if __name__=="__main__":
    main()