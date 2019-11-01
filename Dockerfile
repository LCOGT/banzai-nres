FROM docker.lco.global/banzai:0.27.1-52-g10de024

USER root

COPY --chown=10087:10000 . /lco/banzai-nres

RUN pip install --global-option=build_ext /lco/banzai-nres/ --no-cache-dir

USER archive
