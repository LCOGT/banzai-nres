FROM docker.lco.global/banzai:0.27.5-90-g3fdc726

USER root

RUN conda install -y coveralls sphinx

COPY --chown=10087:10000 . /lco/banzai-nres

RUN pip install --global-option=build_ext /lco/banzai-nres/  halo --no-cache-dir

USER archive
