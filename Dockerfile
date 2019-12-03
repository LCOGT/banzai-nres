FROM docker.lco.global/banzai:0.27.4-67-g47baf06

USER root

RUN conda install -y coveralls sphinx

COPY --chown=10087:10000 . /lco/banzai-nres

RUN pip install --global-option=build_ext /lco/banzai-nres/ --no-cache-dir

USER archive
