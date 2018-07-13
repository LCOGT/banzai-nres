ARG BANZAI_VERSION
FROM docker.lco.global/banzai:${BANZAI_VERSION}

USER root

WORKDIR /lco/banzai-nres

COPY . /lco/banzai-nres

RUN python /lco/banzai-nres/setup.py install

USER archive

ENV HOME /home/archive

WORKDIR /home/archive
