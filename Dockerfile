FROM bentoml/model-server:0.11.0-py312
MAINTAINER ersilia

RUN pip install torch==2.10.0 --index-url https://download.pytorch.org/whl/cpu
RUN pip install torchvision==0.24.1 --index-url https://download.pytorch.org/whl/cpu
RUN pip install transformers==4.51.3
RUN pip install safe-mol==0.1.13

WORKDIR /repo
COPY . /repo
