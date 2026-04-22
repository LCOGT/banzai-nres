FROM ghcr.io/lcogt/banzai:1.35.1

USER root

ENV UV_PROJECT_ENVIRONMENT=/lco/banzai/.venv

COPY pyproject.toml uv.lock /lco/banzai-nres/

RUN uv sync --locked --directory=/lco/banzai-nres --no-install-project

COPY . /lco/banzai-nres

RUN uv sync --locked --directory /lco/banzai-nres

RUN cp /lco/banzai-nres/pytest.ini /home/archive/pytest.ini

RUN chown -R archive:domainusers /home/archive

USER archive
