FROM ghcr.io/lcogt/banzai:1.13.1

USER root

RUN conda install -y coveralls sphinx statsmodels docutils

RUN pip install astropy==5.3.4

COPY --chown=10087:10000 . /lco/banzai-nres

RUN apt-get -y update && apt-get install gcc && \
    pip install /lco/banzai-nres/ --no-cache-dir && \
    apt-get -y remove gcc

RUN chown -R archive /home/archive

USER archive
