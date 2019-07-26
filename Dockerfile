FROM docker.lco.global/banzai:0.26.7-19-g9492cb8

USER root

COPY . /lco/banzai-nres

WORKDIR /lco/banzai-nres

RUN python /lco/banzai-nres/setup.py install

USER archive

ENV HOME /home/archive

WORKDIR /home/archive
