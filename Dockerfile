FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN pip install rdkit==2023.09.3
RUN pip install safe

WORKDIR /repo
COPY . /repo
