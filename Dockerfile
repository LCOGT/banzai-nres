FROM docker.lco.global/banzai:0.26.5-18-ge63e76c

USER root

COPY . /lco/banzai-nres

RUN pip --no-cache-dir install git+https://github.com/lcogt/banzai.git@feature/refactor_context

WORKDIR /lco/banzai-nres

RUN python /lco/banzai-nres/setup.py install

USER archive

ENV HOME /home/archive

WORKDIR /home/archive
