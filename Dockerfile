# Lots of help from the following blog posts:
# https://jcristharif.com/conda-docker-tips.html
# https://opensourcelibs.com/lib/micromamba-docker
# https://uwekorn.com/2021/03/01/deploying-conda-environments-in-docker-how-to-do-it-right.html
FROM mambaorg/micromamba:1.3.0

ENV PYTHONDONTWRITEBYTECODE=true

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml

RUN micromamba install -y --freeze-installed -n base -f /tmp/env.yml just \
    && micromamba clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete

COPY --chown=$MAMBA_USER:$MAMBA_USER . /fastq-dl
WORKDIR /fastq-dl
ENV POETRY_VIRTUALENVS_CREATE=false \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN just install && fastq-dl --version

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]