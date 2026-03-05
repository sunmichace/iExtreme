import os
import sys
import re
import numpy as np
import pandas as pd
import pickle
import uuid
import tempfile
from datetime import datetime
from Bio import SeqIO
from collections import Counter
from sklearn import svm
import subprocess
import warnings
from flask import Flask, request, jsonify, send_from_directory

warnings.filterwarnings("ignore", category=UserWarning, module="sklearn.base")

app = Flask(__name__, static_folder='static', template_folder='templates')
app.config['UPLOAD_FOLDER'] = tempfile.mkdtemp(prefix='iExtreme_')
app.config['MAX_CONTENT_LENGTH'] = 10 * 1024 * 1024
app.config['ALLOWED_EXTENSIONS'] = {'fasta', 'fa', 'faa', 'fna'}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.path.join(SCRIPT_DIR, "optimal")
FEAT_FILE = os.path.join(SCRIPT_DIR, "feat.txt")


def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']


def get_file_size(file_path):
    size = os.path.getsize(file_path)
    if size < 1024 * 1024:
        return f"{size / 1024:.2f}KB"
    return f"{size / (1024 * 1024):.2f}MB"


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


def read_codon_table():
    if not os.path.exists(FEAT_FILE):
        raise Exception(f"特征文件缺失：{FEAT_FILE}")
    try:
        with open(FEAT_FILE, 'r', encoding='utf-8') as f:
            feat = {line.strip(): 0 for line in f if line.strip()}
        return feat
    except Exception as e:
        raise Exception(f"读取特征文件失败：{str(e)}")


# 核心修复：密码子频率计算逻辑（先统计所有序列，再统一计算频率）
def count_freq_codon(seqs):
    codon_table = read_codon_table()
    total_codon = 0  # 统计所有序列的总密码子数
    # 第一步：遍历所有序列，累加密码子计数
    for seq in seqs:
        seq = re.sub(r"\s", "", seq)
        for i in range(0, len(seq), 3):
            codon_2 = seq[i:i + 6]  # 6碱基密码子特征
            if len(codon_2) != 6:
                continue
            if codon_2 in codon_table:
                total_codon += 1
                codon_table[codon_2] += 1
    # 第二步：所有序列统计完成后，统一计算频率（修复缩进）
    if total_codon > 0:
        for k in codon_table:
            codon_table[k] /= total_codon  # 频率 = 该密码子计数 / 总密码子数
    else:
        print("警告：未统计到有效密码子，特征将全部为0")
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

    # 新增：打印特征日志，验证样本差异性（取前3个特征为例）
    feat_sample = list(ccodon_freq_table.items())[:3]
    print(f"样本特征示例（前3个密码子频率）：{dict(feat_sample)}")

    with open(outfile_path, "w") as f:
        feat = list(ccodon_freq_table.keys())
        freq = list(map(str, ccodon_freq_table.values()))
        head_name = "\t".join(feat) + "\tlabel\n"
        freq_line = "\t".join(freq) + "\t-1"
        f.write("refseq\t" + head_name)
        f.write(f"{name}\t{freq_line}\n")
    return outfile_path


class CodonSVM:
    def __init__(self, model_name, test_rate=0.2):
        self.model_name = model_name
        self.model_save_path = os.path.join(MODEL_DIR, f'{self.model_name}.pickle')
        if not os.path.exists(self.model_save_path):
            raise Exception(f"模型文件不存在：{self.model_save_path}")
        self.test_rate = test_rate
        self.model = svm.SVC()

    def load_optimal_model(self):
        try:
            with open(self.model_save_path, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            raise Exception(f"加载模型 {self.model_name} 失败：{str(e)}")

    def predict(self, data_path):
        optimal_model = self.load_optimal_model()
        df = pd.read_csv(data_path, index_col=0, sep="\t")
        X = df.iloc[:, 0:-1].values
        # 数据标准化（使用模型训练时的均值/方差，确保一致性）
        std_x = (X - optimal_model["mean"]) / (optimal_model["var"] **0.5 + 1e-10)
        pred_label = optimal_model["model"].predict(std_x)[0]
        pred_proba = optimal_model["model"].predict_proba(std_x)[0]
        return int(pred_label), np.round(pred_proba, 4)


def run_codon_svm(data_path, model_name):
    codon_svm = CodonSVM(model_name)
    return codon_svm.predict(data_path)


@app.route('/')
def index():
    return send_from_directory(app.template_folder, 'index.html')


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

        models = [{"name": "exheat"}, {"name": "exsalt"}, {"name": "expH"}]
        task_id = f"task_{datetime.now().strftime('%Y%m%d%H%M%S')}_{uuid.uuid4().hex[:8]}"
        results = {
            "status": "success",
            "task_id": task_id,
            "file_results": [],
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
                "models": {}
            }

            try:
                # 步骤1：Prodigal预测
                prodigal_out = os.path.join(result_dir, f"prodigal_{os.path.splitext(file_name)[0]}")
                cds_file, _ = run_prodigal(file_path, prodigal_out)

                # 步骤2：预处理（新增特征日志）
                print(f"\n=== 文件 {file_name} 特征提取 ===")
                process_file = os.path.join(result_dir, f"process_{os.path.splitext(file_name)[0]}.tsv")
                preprocess_data(cds_file, process_file)
                seqs = get_one_seqs(cds_file)
                file_result["sequence_count"] = len(seqs) if seqs[0] != "error\n" else 0
                print(f"文件 {file_name} 包含序列数：{file_result['sequence_count']}")

                # 步骤3：模型预测
                print("\n" + "="*80)
                print(f"文件 {file_name} 的三个模型预测结果：")
                print("="*80)

                for model in models:
                    model_name = model["name"]
                    print(f"\n【开始运行模型：{model_name}】")

                    try:
                        pred_label, pred_proba = run_codon_svm(process_file, model_name)
                        n_classes = len(pred_proba)
                        model_result = {"label": int(pred_label), "probabilities": {}}
                        for i in range(n_classes):
                            model_result["probabilities"][f"proba_label_{i}"] = float(pred_proba[i])

                        print(f"【模型 {model_name} 运行成功】")
                        print(f"预测标签：{pred_label}")
                        print(f"概率分布（共{ n_classes }类）：")
                        for i in range(n_classes):
                            print(f"  proba_label_{i}: {pred_proba[i]:.4f}")

                        file_result["models"][model_name] = model_result

                    except Exception as e:
                        error_msg = f"模型 {model_name} 运行失败：{str(e)}"
                        file_result["models"][model_name] = {"error": error_msg}
                        results["warnings"].append(error_msg)
                        print(f"【模型 {model_name} 运行失败】：{error_msg}")
                        continue

                print("\n" + "="*80 + "\n")
                results["file_results"].append(file_result)

            except Exception as e:
                file_result["status"] = "error"
                file_result["error"] = str(e)
                results["file_results"].append(file_result)
                results["warnings"].append(f"文件 {file_name} 预处理失败：{str(e)}")
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
    os.makedirs(os.path.join(SCRIPT_DIR, "data"), exist_ok=True)
    app.run(
        host='0.0.0.0',
        port=9060,
        debug=False,
        threaded=True
    )