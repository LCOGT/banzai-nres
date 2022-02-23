FROM docker.lco.global/banzai:1.9.5

USER root

RUN conda install -y coveralls sphinx statsmodels docutils=0.15

RUN pip install astropy==4.2

COPY --chown=10087:10000 . /lco/banzai-nres

RUN apt-get -y update && apt-get install gcc && \
    pip install /lco/banzai-nres/ --no-cache-dir && \
    apt-get -y remove gcc

RUN chown -R archive /home/archive

USER archive
