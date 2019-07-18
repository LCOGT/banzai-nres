FROM docker.lco.global/banzai:0.26.3
USER root

RUN pip install --no-cache-dir  git+https://github.com/lcogt/banzai.git@fix/circular_import

WORKDIR /lco/banzai-nres

COPY . /lco/banzai-nres

RUN python /lco/banzai-nres/setup.py install

USER archive

ENV HOME /home/archive

WORKDIR /home/archive
