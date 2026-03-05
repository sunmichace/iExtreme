iExtreme - 极端微生物生存条件预测工具
项目概述
iExtreme 是一款基于机器学习（CodonSVM + TSpH 模型）的极端微生物生存条件预测工具，支持上传 FASTA 格式文件，自动完成基因预测、特征提取、模型推理，最终输出微生物的最适温度、最适盐度、最适 pH 值，以及热 / 盐 /pH 耐受性模型的预测结果（其中温度模型返回 label1+label2 的概率和）。
核心功能
支持上传 FASTA/fa/faa/fna 格式文件（最大 10MB）
基于 Prodigal 完成基因预测（CDS / 蛋白质序列生成）
CodonSVM 模型预测：exheat（热耐受）、exsalt（盐耐受）、expH（pH 耐受）
特殊处理：exheat 模型返回 label1+label2 的概率和（非固定 label2）
exsalt/expH 模型保持固定 label=1
TSpH 模型预测最适温度、盐度、pH 值，及对应的差异分析数据
自动清理临时文件，保证服务器存储整洁
环境依赖
1. 必备软件
Prodigal（基因预测工具）：用于从 FASTA 文件提取 CDS 和蛋白质序列
bash
运行
# 通过conda安装（推荐）
conda install -c bioconda prodigal
2. Python 依赖包
创建requirements.txt文件，内容如下：
txt
pandas>=1.5.0
numpy>=1.21.0
scikit-learn>=1.0.0
biopython>=1.79
flask>=2.0.0
joblib>=1.1.0
multiprocessing>=2.6.0.1
uuid>=1.30
tempfile>=0.1.0
re>=2.2.1
subprocess32>=3.5.4  # Windows系统可省略
安装依赖：
bash
运行
pip install -r requirements.txt
项目结构
plaintext
iExtreme/
├── app_all.py               # 主程序（Flask后端）
├── static/              # 静态资源目录（存放echarts.min.js等）
├── templates/           # 前端页面目录（存放all.html）
├── optimal/             # 模型目录（存放CodonSVM和TSpH模型文件）
│   ├── exheat.pickle    # 热耐受模型
│   ├── exsalt.pickle    # 盐耐受模型
│   ├── expH.pickle      # pH耐受模型
│   ├── train_temp.pkl   # TSpH温度预测模型
│   ├── train_salt.pkl   # TSpH盐度预测模型
│   ├── train_pH.pkl     # TSpH pH预测模型
│   └── 对应模型的.f特征文件
├── data/                # 差异数据目录
│   ├── my_pred_train_df_temp.csv.diff  # 温度差异数据
│   ├── my_pred_train_df_salt.csv.diff  # 盐度差异数据
│   └── my_pred_train_df_pH.csv.diff    # pH差异数据
├── feat.txt             # CodonSVM特征文件（6联密码子特征）
└── README.md            # 项目说明文档
快速开始
1. 准备工作
将训练好的模型文件（.pickle/.pkl/.f）放入optimal目录
将差异数据文件放入data目录
确保feat.txt特征文件存在且格式正确（每行一个 6 联密码子特征）
2. 启动服务
# 运行主程序
python app_all.py
服务启动后，默认监听：http://0.0.0.0:9070
3. 使用流程
打开浏览器访问 http://localhost:9070
点击「Select File」上传 FASTA 格式文件（最大 10MB）
点击「Start Prediction」开始预测
等待预测完成后，查看：
最适温度 / 盐度 /pH 值
三个模型的预测分数（exheat 为 label1+2 概率和）
温度 / 盐度 /pH 的差异分析图表
API 接口说明
1. 文件上传接口
路径：/upload
请求方式：POST
请求体：form-data 格式，包含file字段（上传的文件）
返回格式（JSON）：
json
{
  "status": "success",  // success/error
  "filePath": "临时文件路径",
  "filename": "原始文件名",
  "size": "文件大小（如1.25MB）",
  "msg": "错误信息（仅error时返回）"
}
2. 预测处理接口
路径：/process_files
请求方式：POST
请求头：Content-Type: application/json
请求体（JSON）：
json
{
  "filePaths": ["上传返回的临时文件路径"]
}
返回格式（JSON）：
json
{
  "status": "success",
  "task_id": "任务ID",
  "file_results": [
    {
      "filename": "文件名",
      "size": "文件大小",
      "status": "success",
      "sequence_count": 序列数量,
      "models": {
        "exheat": {
          "label": 0.85, 
          "probabilities": {
            "proba_label_0": 0.1,
            "proba_label_1": 0.4,
            "proba_label_2": 0.45
          }
        },
        "exsalt": {
          "label": 1,
          "probabilities": {...}
        },
        "expH": {
          "label": 1,
          "probabilities": {...}
        }
      },
      "temperature_opt": 55.2,  // 最适温度
      "salinity_opt": 3.5,      // 最适盐度
      "ph_opt": 7.2,            // 最适pH
      "temp_x": "差异数据x",
      "temp_y": "差异数据y",
      "salinity_x": "差异数据x",
      "salinity_y": "差异数据y",
      "ph_x": "差异数据x",
      "ph_y": "差异数据y"
    }
  ],
  "warnings": []  // 警告信息列表
}
关键配置说明
表格
配置项	说明	默认值
MAX_CONTENT_LENGTH	上传文件大小限制	10MB
ALLOWED_EXTENSIONS	支持的文件格式	fasta/fa/faa/fna
MODEL_DIR	模型文件目录	./optimal
DATA_DIR	差异数据目录	./data
UPLOAD_FOLDER	临时文件存储目录（自动生成）	tempfile.mkdtemp(prefix='iExtreme_')
服务端口	Flask 运行端口	9070
注意事项
模型文件要求：
CodonSVM 模型为 pickle 格式，需包含mean/var/model字段
TSpH 模型为 joblib 格式，需配套.f特征文件（包含 mean/std/features）
临时文件：预测完成后会自动删除临时目录和上传文件，无需手动清理
编码问题：特征文件 / CSV 文件建议使用utf-8-sig编码，避免中文 / BOM 头问题
并发处理：服务启用threaded=True支持多线程，可根据服务器配置调整
兼容性：Prodigal 在 Windows 系统下可能需要手动配置环境变量
常见问题
Q1: 提示「Prodigal 未安装」
A1: 确保通过 conda 安装 Prodigal，且已加入系统环境变量，可在终端执行prodigal -h验证。
Q2: 温度模型仍返回 label2
A2: 检查process_files函数中exheat模型的处理逻辑，确保：
已删除adjusted_label = 2硬编码
正确计算sum_label_1_2 = pred_proba[1] + pred_proba[2]
前端已修改为读取exheat.label而非proba_label_2
Q3: 模型加载失败
A3: 检查模型文件路径是否正确，模型文件是否完整，pickle/joblib 版本是否兼容。
许可证
本项目仅供学术研究使用，如需商用请联系开发者。