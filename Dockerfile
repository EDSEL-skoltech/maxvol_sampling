FROM python:3.6

LABEL maintainer='Mikhail.Gasanov@skoltech.ru'

RUN DEBIAN_FRONTEND='noninteractive' apt-get update && apt-get upgrade -y --no-install-recommends && \
    DEBIAN_FRONTEND='noninteractive' apt-get install -y --no-install-recommends liblapack-dev liblapack3 libopenblas-base

COPY requirements.txt .

RUN pip install -r requirements.txt && pip install maxvolpy
