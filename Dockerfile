FROM julia:1.12-bookworm

ENV JULIA_DEPOT_PATH=/usr/local/julia-depot
ENV JULIA_PKG_PRECOMPILE_AUTO=0

RUN apt-get update \
    && apt-get install -y --no-install-recommends git build-essential gcc \
    && rm -rf /var/lib/apt/lists/*

RUN useradd --create-home --shell /bin/bash --uid 1000 vscode

WORKDIR /home/vscode

COPY start-pluto.sh ./

WORKDIR /tmp/julia-project

COPY Project.toml Manifest.toml ./

# 1. Installation globale des paquets
RUN mkdir -p "${JULIA_DEPOT_PATH}" \
    && julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using Plots, AlgebraicSolving, Pluto, IJulia'

# 2. Création de l'image système avec les paquets pré-compilés
RUN julia -e 'using PackageCompiler; create_sysimage(["Plots", "AlgebraicSolving", "Pluto", "IJulia"], sysimage_path="/usr/local/julia-depot/sys_albatross.so")'

# 3. Création d'un noyau Jupyter spécifique utilisant cette image
RUN julia -e 'using IJulia; installkernel("Julia ALBATROSS", "--sysimage=/usr/local/julia-depot/sys_albatross.so"); notebook()'

RUN chmod -R a+rwX "${JULIA_DEPOT_PATH}"
