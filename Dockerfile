FROM julia:1.12-bookworm

ENV JULIA_DEPOT_PATH=/usr/local/julia-depot
ENV JULIA_PKG_PRECOMPILE_AUTO=0

# System dependencies + Python/Jupyter
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git \
        build-essential \
        python3 \
        python3-pip \
        python3-venv \
	htop \
    && rm -rf /var/lib/apt/lists/*

# Install Jupyter independently of Julia/Conda
RUN python3 -m venv /opt/jupyter \
    && /opt/jupyter/bin/pip install --no-cache-dir \
        jupyter \
        notebook

ENV PATH="/opt/jupyter/bin:${PATH}"
ENV JUPYTER="/opt/jupyter/bin/jupyter"

# Create user
RUN useradd --create-home --shell /bin/bash --uid 1000 vscode

WORKDIR /tmp/julia-project

COPY Project.toml Manifest.toml ./

RUN mkdir -p "${JULIA_DEPOT_PATH}" \
    && julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using Plots, AlgebraicSolving, Pluto' \
    && chmod -R a+rwX "${JULIA_DEPOT_PATH}"

USER julia

WORKDIR /workspace

# Install IJulia and explicitly tell it which Jupyter to use.
RUN julia --project=. -e 'ENV["JUPYTER"] = "/opt/jupyter/bin/jupyter"; using Pkg; Pkg.build("IJulia")'
CMD ["julia", "-e", "using IJulia; notebook()"]
