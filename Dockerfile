FROM continuumio/miniconda3:latest
WORKDIR /app
COPY environment.yml .
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ && \
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
conda config --set show_channel_urls yes && \
conda env create -f environment.yml && \
conda clean -a -y && \
conda run -n fastani_prodigal pip install --upgrade pip
COPY . .
EXPOSE 9060
CMD ["conda", "run", "--no-capture-output", "-n", "fastani_prodigal", "python", "app_salt.py"]
