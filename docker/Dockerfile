FROM ubuntu:15.10

# File author / Maintainer
MAINTAINER Charles Markello <cmarkell@ucsc.edu>

RUN apt-get update && apt-get install -y \
    git \
    python-dev \
    python-pip \
    wget \
    curl \
    apt-transport-https \
    ca-certificates

# Get the Docker binary
RUN curl https://get.docker.com/builds/Linux/x86_64/docker-1.12.3.tgz \
    | tar -xvzf - --transform='s,[^/]*/,,g' -C /usr/local/bin/ \
    && chmod u+x /usr/local/bin/docker

# Install Toil
RUN pip install toil==3.6.0

# Install toil-vg
RUN pip install toil-vg==1.1.0a1.dev312

COPY wrapper.sh /opt/vg-pipeline/

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["sh", "/opt/vg-pipeline/wrapper.sh"]
CMD ["--help"]
