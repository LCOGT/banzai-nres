FROM docker.lco.global/banzai:0.28.9-196-g370048b

USER root

RUN conda install -y coveralls sphinx statsmodels docutils=0.15

COPY --chown=10087:10000 . /lco/banzai-nres

RUN pip install --global-option=build_ext /lco/banzai-nres/ --no-cache-dir

USER archive
