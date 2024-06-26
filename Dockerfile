FROM bentoml/model-server:0.11.0-py39
MAINTAINER ersilia


RUN pip install protobuf==3.18.3
RUN pip install tqdm==4.66.1
RUN pip install accelerate==0.25.0
RUN pip install datamol==0.12.2
RUN pip install datasets==2.16.1
RUN pip install evaluate==0.4.1
RUN pip install tokenizers==0.15.0
RUN pip install transformers==4.36.2
RUN pip install typer==0.9.0
RUN pip install universal-pathlib==0.1.4
RUN pip install wandb==0.16.0
RUN pip install rdkit==2023.9.3
RUN pip install pandas==2.1
RUN pip install --no-deps safe-mol==0.1.4
RUN conda install -c conda-forge xorg-libxrender xorg-libxtst

WORKDIR /repo
COPY . /repo
