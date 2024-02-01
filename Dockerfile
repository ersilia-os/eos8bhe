FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN pip install safe-mol==0.1.4

WORKDIR /repo
COPY . /repo
