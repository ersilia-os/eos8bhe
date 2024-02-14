FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN pip install --no-deps safe-mol==0.1.4
RUN pip install tqdm==4.66.1
RUN pip install accelerate==0.25.0
RUN pip install datamol==0.12.2
RUN pip install datasets==2.16.1
RUN pip install evaluate==0.4.1
RUN pip install tokenizers==0.15.0
RUN pip install transformers==4.36.2
RUN pip install typer==0.9.0
RUN pip install universal-pathlib==0.1.4

WORKDIR /repo
COPY . /repo
