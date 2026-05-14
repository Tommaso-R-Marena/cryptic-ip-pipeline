FROM mambaorg/micromamba:1.5.8

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl git build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/crypticip
COPY environment.yml /tmp/environment.yml
RUN micromamba env create -y -n crypticip -f /tmp/environment.yml \
    && micromamba clean -a -y
ENV PATH=/opt/conda/envs/crypticip/bin:$PATH

COPY . /opt/crypticip
RUN /opt/conda/envs/crypticip/bin/pip install --no-deps -e .

ENV PYTHONUNBUFFERED=1
CMD ["crypticip", "--help"]
