FROM docker.lco.global/banzai:0.26.5-8-g15d6b7d

USER root

WORKDIR /lco/banzai-nres

COPY . /lco/banzai-nres

RUN pip install --no-cache-dir --upgrade git+https://github.com/lcogt/banzai.git@fix/circular_import

RUN python /lco/banzai-nres/setup.py install

USER archive

ENV HOME /home/archive

WORKDIR /home/archive
