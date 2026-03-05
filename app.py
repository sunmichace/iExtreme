import os
import sys
import re
import joblib
import numpy as np
import pandas as pd
import pickle
import uuid
import tempfile
from datetime import datetime
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from collections import Counter
from sklearn import svm
from sklearn.model_selection import GridSearchCV
import subprocess
from flask import Flask, request, jsonify, send_from_directory

# 初始化 Flask 应用
app = Flask(__name__, static_folder='static', template_folder='templates')
app.config['UPLOAD_FOLDER'] = tempfile.mkdtemp(prefix='iExtreme_')
app.config['MAX_CONTENT_LENGTH'] = 10 * 1024 * 1024  # 10MB限制
app.config['ALLOWED_EXTENSIONS'] = {'fasta', 'fa', 'faa', 'fna'}

# 路径配置
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.path.join(SCRIPT_DIR, "optimal")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
FEAT_FILE = os.path.join(SCRIPT_DIR, "feat.txt")
SEED_FILE = os.path.join(SCRIPT_DIR, "seed_pool.txt")


# 基础工具函数（不变）
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']


def get_file_size(file_path):
    size = os.path.getsize(file_path)
    if size < 1024 * 1024:
        return f"{size / 1024:.2f}KB"
    return f"{size / (1024 * 1024):.2f}MB"


# 1. Prodigal 基因预测（不变）
def run_prodigal(input_fasta, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    cds_output = os.path.join(output_dir, "user_pred.cds")
    pep_output = os.path.join(output_dir, "user_pred.pep")
    gbk_output = os.path.join(output_dir, "prodigal.gbk")

    cmd = [
        "prodigal",
        "-i", input_fasta,
        "-o", gbk_output,
        "-a", pep_output,
        "-d", cds_output,
        "-p", "meta"
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return cds_output, pep_output
    except subprocess.CalledProcessError as e:
        raise Exception(f"Prodigal运行失败：{e.stderr}")
    except FileNotFoundError:
        raise Exception("Prodigal未安装！请通过 'conda install -c bioconda prodigal' 安装")


# 2. 数据预处理（不变）
def read_codon_table():
    if not os.path.exists(FEAT_FILE):
        raise Exception(f"特征文件缺失：{FEAT_FILE}")
    try:
        with open(FEAT_FILE, 'r', encoding='utf-8') as f:
            feat = {line.strip(): 0 for line in f if line.strip()}
        if len(feat) == 0:
            raise Exception("特征文件为空")
        return feat
    except Exception as e:
        raise Exception(f"读取特征文件失败：{str(e)}")


def count_freq_codon(seqs):
    codon_table = read_codon_table()
    total_codon = 0
    for seq in seqs:
        seq = re.sub(r"\s", "", seq)
        for i in range(0, len(seq), 3):
            codon_2 = seq[i:i + 6]
            if len(codon_2) != 6:
                continue
            if codon_2 in codon_table:
                total_codon += 1
                codon_table[codon_2] += 1
    if total_codon > 0:
        for k in codon_table:
            codon_table[k] /= total_codon
    return codon_table


def get_one_seqs(path):
    try:
        with open(path, 'r') as f:
            c = f.read()
        return re.findall(">.+\n([^>]+)", c)
    except:
        return ["error\n"]


def preprocess_data(cds_file, outfile_path):
    name = "user"
    seqs = get_one_seqs(cds_file)
    if seqs and seqs[0] == "error\n":
        ccodon_freq_table = {"error": "error"}
    else:
        ccodon_freq_table = count_freq_codon(seqs)

    with open(outfile_path, "w") as f:
        feat = list(ccodon_freq_table.keys())
        freq = list(map(str, ccodon_freq_table.values()))
        head_name = "\t".join(feat) + "\tlabel\n"
        freq_line = "\t".join(freq) + "\t-1"
        f.write("refseq\t" + head_name)
        f.write(f"{name}\t{freq_line}\n")
    return outfile_path


# 3. CodonSVM 预测（保留文件级预测）
class CodonSVM:
    def __init__(self, model_name, optimal_label):
        self.model_name = model_name
        self.model = svm.SVC(probability=True)
        self.optimal_label = optimal_label
        self.model_save_path = os.path.join(MODEL_DIR, f'{self.model_name}.pickle')
        self.seed_pool = self.load_seed_pool()

    def load_seed_pool(self):
        if not os.path.exists(SEED_FILE):
            raise Exception(f"种子池文件缺失：{SEED_FILE}")
        try:
            with open(SEED_FILE, "r") as f:
                return [int(seed.strip()) for seed in f if seed.strip()]
        except Exception as e:
            raise Exception(f"读取种子池失败：{str(e)}")

    def load_optimal_model(self):
        if not os.path.exists(self.model_save_path):
            raise Exception(f"SVM模型缺失：{self.model_save_path}")
        try:
            with open(self.model_save_path, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            raise Exception(f"加载SVM模型失败：{str(e)}")

    def predict(self, data_path):
        optimal_model = self.load_optimal_model()
        df = pd.read_csv(data_path, index_col=0, sep="\t")
        X = df.iloc[:, 0:-1].values
        std_x = (X - optimal_model["mean"]) / (optimal_model["var"] ** 0.5 + 1e-10)
        pred_label = optimal_model["model"].predict(std_x)[0]
        pred_proba = optimal_model["model"].predict_proba(std_x)[0]
        return int(pred_label), np.round(pred_proba, 4)


def run_codon_svm(data_path, model_name, optimal_label):
    codon_svm = CodonSVM(model_name, optimal_label)
    return codon_svm.predict(data_path)


# 4. TSpH 数值预测（完整返回差异数据）
def load_means_stds(predictor):
    means = dict()
    stds = dict()
    features = list()
    if not os.path.exists(predictor):
        raise Exception(f"TSpH特征文件缺失：{predictor}")
    with open(predictor, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cont = line.strip().split()
            if len(cont) != 3:
                continue
            means[cont[0]] = float(cont[1])
            stds[cont[0]] = float(cont[2])
            features.append(cont[0])
    return means, stds, features


def load_model(path):
    if not os.path.exists(path):
        raise Exception(f"TSpH模型缺失：{path}")
    try:
        model = joblib.load(path)
        feature_path = path.replace('.pkl', '.f')
        means, stds, features = load_means_stds(feature_path)
        return model, means, stds, features
    except Exception as e:
        raise Exception(f"加载TSpH模型失败：{str(e)}")


def do_count(seq):
    dimers = Counter()
    for i in range(len(seq) - 1):
        dimer = seq[i:i + 2]
        dimers[dimer] += 1.0
    return dimers


def count_dimer(fasta_file, p=2):
    seqs = [str(rec.seq).upper() for rec in SeqIO.parse(fasta_file, 'fasta')]
    if not seqs:
        raise Exception("蛋白质文件中无有效序列")
    num_cpus = min(int(p), cpu_count(), len(seqs))
    with Pool(num_cpus) as pool:
        results = pool.map(do_count, seqs)
    return sum(results, Counter()), seqs


def get_dimer_frequency(fasta_file, p=2):
    total_dimers, seqs = count_dimer(fasta_file, p)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
    dimers_fq = {a1 + a2: total_dimers.get(a1 + a2, 0.0) for a1 in amino_acids for a2 in amino_acids}
    total = sum(dimers_fq.values())
    if total == 0:
        raise Exception("未检测到有效蛋白质二聚体")
    for key in dimers_fq:
        dimers_fq[key] /= total
    return dimers_fq, seqs


def predict_tsph(fasta_file, model, means, stds, features, p=2):
    dimers_fq, seqs = get_dimer_frequency(fasta_file, p)
    Xs = [(dimers_fq[fea] - means[fea]) / stds[fea] for fea in features]
    Xs = np.array(Xs).reshape([1, len(Xs)])
    pred = model.predict(Xs)[0]
    return np.round(pred, 2), seqs


def get_diff_x_y(y, diff_path):
    if not os.path.exists(diff_path):
        raise Exception(f"差异文件缺失：{diff_path}")
    try:
        df = pd.read_csv(diff_path)
        diff_x = df["diff_x"].values + y
        diff_x[diff_x < 0] = 0
        diff_y = df["diff_p"].values
        return np.round(diff_x, 2), np.round(diff_y, 2)
    except Exception as e:
        raise Exception(f"读取差异文件失败：{str(e)}")


def run_predict_tsph(pep_file):
    model_paths = {
        "Temp": os.path.join(MODEL_DIR, "train_temp.pkl"),
        "Salinity": os.path.join(MODEL_DIR, "train_salt.pkl"),
        "pH": os.path.join(MODEL_DIR, "train_pH.pkl")
    }
    diff_paths = {
        "Temp": os.path.join(DATA_DIR, "my_pred_train_df_temp.csv.diff"),
        "Salinity": os.path.join(DATA_DIR, "my_pred_train_df_salt.csv.diff"),
        "pH": os.path.join(DATA_DIR, "my_pred_train_df_pH.csv.diff")
    }

    temp_model, means_T, stds_T, features_T = load_model(model_paths["Temp"])
    salinity_model, means_S, stds_S, features_S = load_model(model_paths["Salinity"])
    ph_model, means_pH, stds_pH, features_pH = load_model(model_paths["pH"])

    temp_pred, seqs = predict_tsph(pep_file, temp_model, means_T, stds_T, features_T)
    salinity_pred, _ = predict_tsph(pep_file, salinity_model, means_S, stds_S, features_S)
    ph_pred, _ = predict_tsph(pep_file, ph_model, means_pH, stds_pH, features_pH)

    temp_diff_x, temp_diff_y = get_diff_x_y(temp_pred, diff_paths["Temp"])
    salinity_diff_x, salinity_diff_y = get_diff_x_y(salinity_pred, diff_paths["Salinity"])
    ph_diff_x, ph_diff_y = get_diff_x_y(ph_pred, diff_paths["pH"])

    # 格式化差异数据为字符串（与原始输出一致）
    return {
        "temp": {
            "pred": temp_pred,
            "diff_x": ",".join(map(str, temp_diff_x)),
            "diff_y": ",".join(map(str, temp_diff_y))
        },
        "salinity": {
            "pred": salinity_pred,
            "diff_x": ",".join(map(str, salinity_diff_x)),
            "diff_y": ",".join(map(str, salinity_diff_y))
        },
        "ph": {
            "pred": ph_pred,
            "diff_x": ",".join(map(str, ph_diff_x)),
            "diff_y": ",".join(map(str, ph_diff_y))
        },
        "seq_count": len(seqs)
    }


# Flask API 接口（核心修改：只返回文件级结果）
@app.route('/')
def index():
    return send_from_directory(app.template_folder, 'all.html')


@app.route('/upload', methods=['POST'])
def upload_file():
    try:
        if 'file' not in request.files:
            return jsonify({"status": "error", "msg": "未检测到文件"})
        file = request.files['file']
        if file.filename == '':
            return jsonify({"status": "error", "msg": "未选择文件"})
        if not allowed_file(file.filename):
            return jsonify({"status": "error", "msg": f"不支持的格式：{file.filename}"})

        file_id = f"file_{datetime.now().strftime('%Y%m%d%H%M%S')}_{uuid.uuid4().hex[:8]}"
        filename = f"{file_id}_{file.filename}"
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(file_path)

        return jsonify({
            "status": "success",
            "filePath": file_path,
            "filename": file.filename,
            "size": get_file_size(file_path)
        })
    except Exception as e:
        return jsonify({"status": "error", "msg": f"文件上传失败：{str(e)}"})


@app.route('/process_files', methods=['POST'])
def process_files():
    try:
        data = request.get_json()
        file_paths = data.get('filePaths')
        if not file_paths or not isinstance(file_paths, list):
            return jsonify({"status": "error", "msg": "请传入有效的文件路径列表"})

        task_id = f"task_{datetime.now().strftime('%Y%m%d%H%M%S')}_{uuid.uuid4().hex[:8]}"
        results = {
            "status": "success",
            "task_id": task_id,
            "file_results": [],  # 存储每个文件的整体预测结果
            "warnings": []
        }

        result_dir = os.path.join(app.config['UPLOAD_FOLDER'], task_id)
        os.makedirs(result_dir, exist_ok=True)

        for file_path in file_paths:
            file_name = os.path.basename(file_path)
            file_size = get_file_size(file_path)
            file_result = {
                "filename": file_name,
                "size": file_size,
                "status": "success",
                "sequence_count": 0,
                "error": "",
                # 预测结果字段
                "extreme_heat": "No",
                "extreme_heat_proba": 0.0,
                "extreme_salt": "No",
                "extreme_salt_proba": 0.0,
                "extreme_ph": "No",
                "extreme_ph_proba": 0.0,
                "temperature_opt": 0.0,
                "salinity_opt": 0.0,
                "ph_opt": 0.0,
                # 详细差异数据
                "temp_x": "",
                "temp_y": "",
                "salinity_x": "",
                "salinity_y": "",
                "ph_x": "",
                "ph_y": ""
            }

            try:
                # 步骤1：Prodigal预测
                prodigal_out = os.path.join(result_dir, f"prodigal_{os.path.splitext(file_name)[0]}")
                cds_file, pep_file = run_prodigal(file_path, prodigal_out)

                # 步骤2：预处理
                process_file = os.path.join(result_dir, f"process_{os.path.splitext(file_name)[0]}.tsv")
                preprocess_data(cds_file, process_file)
                seqs = get_one_seqs(cds_file)
                file_result["sequence_count"] = len(seqs) if seqs[0] != "error\n" else 0

                # 步骤3：SVM预测（文件级）
                heat_label, heat_proba = run_codon_svm(process_file, "exheat", "2")
                salt_label, salt_proba = run_codon_svm(process_file, "exsalt", "2")
                ph_label, ph_proba = run_codon_svm(process_file, "expH", "2")

                file_result["extreme_heat_proba"] = float(heat_proba[-1])  # 取最后一个元素
                file_result["extreme_salt_proba"] = float(salt_proba[-1])
                file_result["extreme_ph_proba"] = float(ph_proba[-1])

                # 步骤4：TSpH预测（包含详细差异数据）
                tsph_result = run_predict_tsph(pep_file)
                file_result["temperature_opt"] = float(tsph_result["temp"]["pred"])
                file_result["salinity_opt"] = float(tsph_result["salinity"]["pred"])
                file_result["ph_opt"] = float(tsph_result["ph"]["pred"])
                # 详细差异数据（与原始输出一致）
                file_result["temp_x"] = tsph_result["temp"]["diff_x"]
                file_result["temp_y"] = tsph_result["temp"]["diff_y"]
                file_result["salinity_x"] = tsph_result["salinity"]["diff_x"]
                file_result["salinity_y"] = tsph_result["salinity"]["diff_y"]
                file_result["ph_x"] = tsph_result["ph"]["diff_x"]
                file_result["ph_y"] = tsph_result["ph"]["diff_y"]

                # 打印原始格式输出（后端控制台）
                print(
                    f"原始预测结果 - 温度：{file_result['temperature_opt']}，盐度：{file_result['salinity_opt']}，pH：{file_result['ph_opt']}")
                print(f">Temp_x={file_result['temp_x']}")
                print(f">Temp_y={file_result['temp_y']}")
                print(f">Salinity_x={file_result['salinity_x']}")
                print(f">Salinity_y={file_result['salinity_y']}")
                print(f">pH_x={file_result['ph_x']}")
                print(f">pH_y={file_result['ph_y']}")

                results["file_results"].append(file_result)

            except Exception as e:
                file_result["status"] = "error"
                file_result["error"] = str(e)
                results["file_results"].append(file_result)
                results["warnings"].append(f"文件 {file_name} 处理失败：{str(e)}")
                continue

        # 清理临时文件
        import shutil
        shutil.rmtree(result_dir)
        for file_path in file_paths:
            os.remove(file_path)

        return jsonify(results)

    except Exception as e:
        return jsonify({"status": "error", "msg": f"预测流程失败：{str(e)}"})


if __name__ == "__main__":
    os.makedirs(MODEL_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)
    app.run(
        host='0.0.0.0',
        port=9070,
        debug=False,
        threaded=True
    )