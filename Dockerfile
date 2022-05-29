# Lots of help from the following blog posts:
# https://jcristharif.com/conda-docker-tips.html
# https://opensourcelibs.com/lib/micromamba-docker
# https://uwekorn.com/2021/03/01/deploying-conda-environments-in-docker-how-to-do-it-right.html
FROM mambaorg/micromamba:0.23.3

ENV PYTHONDONTWRITEBYTECODE=true

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml

RUN micromamba install -y --freeze-installed -n base -f /tmp/env.yml \
    && micromamba clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete

COPY --chown=$MAMBA_USER:$MAMBA_USER fastq-dl.py /bin/fastq-dl
