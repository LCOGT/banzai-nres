FROM docker.lco.global/banzai:0.27.4-59-g53cca61

USER root

RUN conda install -y coveralls sphinx

COPY --chown=10087:10000 . /lco/banzai-nres

RUN pip install --global-option=build_ext /lco/banzai-nres/ --no-cache-dir

USER archive
