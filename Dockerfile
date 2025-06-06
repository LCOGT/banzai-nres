FROM ghcr.io/lcogt/banzai:1.24.1

USER root

RUN poetry config virtualenvs.create false

COPY pyproject.toml poetry.lock /lco/banzai-nres/

RUN poetry install --directory=/lco/banzai-nres -E cpu --no-root --no-cache

COPY . /lco/banzai-nres

RUN poetry install --directory /lco/banzai-nres -E cpu --no-cache

RUN cp /lco/banzai-nres/pytest.ini /home/archive/pytest.ini

RUN chown -R archive:domainusers /home/archive

USER archive
