FROM python:3.6.6-slim-jessie
ENTRYPOINT  ["/bin/bash", "-c", "while true; do sleep 100; done"]
USER root

RUN apt-get update && apt-get -y --no-install-recommends install build-essential git && rm -rf /var/lib/apt/lists/*

RUN pip install numpy Cython && rm -rf ~/.cache/pip

RUN git clone https://github.com/kbarbary/sep.git /usr/src/sep
WORKDIR /usr/src/sep
RUN python setup.py install

WORKDIR /lco/banzai-nres

COPY . /lco/banzai-nres

RUN python /lco/banzai-nres/setup.py install

USER archive

ENV HOME /home/archive

WORKDIR /home/archive
