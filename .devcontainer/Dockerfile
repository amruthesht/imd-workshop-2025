FROM ghcr.io/becksteinlab/streaming-md-docker:main-common-cpu

COPY env.yaml env.yaml
RUN . /opt/conda/etc/profile.d/conda.sh && \
        conda env create -n workshop --file=env.yaml

RUN  . /opt/conda/etc/profile.d/conda.sh && \
        conda clean --all -y -f

# Install additional tools
RUN apt-get update && apt-get install -y \
    zsh \
    vim \
    less \
    curl \
    wget \
    git \
    tree \
    htop \
    neofetch \
    bat \
    fd-find \
    ripgrep \
    fonts-powerline \
    && rm -rf /var/lib/apt/lists/*

COPY create_conda_user /opt/docker/bin/create_conda_user
RUN chmod +x /opt/docker/bin/create_conda_user

# Set GITHUB_USER environment variable for container
ARG GITHUB_USER
ENV GITHUB_USER=${GITHUB_USER}

RUN /opt/docker/bin/create_conda_user

# Codespaces does not use these!
# ENTRYPOINT [ "/opt/conda/bin/tini", "--", "/opt/docker/bin/workshop_entrypoint" ]
# CMD [ "/bin/zsh" ]

USER ${GITHUB_USER:-workshop-user} 
