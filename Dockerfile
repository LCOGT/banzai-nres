FROM docker.lco.global/banzai:1.0.5

USER root

RUN conda install -y coveralls sphinx statsmodels docutils=0.15

RUN pip install astropy==4.1rc1

COPY --chown=10087:10000 . /lco/banzai-nres

RUN pip install /lco/banzai-nres/ --no-cache-dir

RUN chown -R archive /home/archive

USER archive
