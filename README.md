# ALBATROSS 2026 — Practical Session

This repository contains the material for the practical session of the **ALBATROSS 2026 workshop** (Algebraic Methods for Polynomial System Solving).

The practical session uses [Julia](https://julialang.org/) and [AlgebraicSolving.jl](https://algebraic-solving.github.io/).

## Repository contents

* Jupyter notebooks (`*.ipynb`) for the practical session.
* `solutions/`: solutions for the practical session.
* `julia_files/`: Julia versions of the notebooks.
* `install.sh`: installation script for a local Julia installation.
* `install_docker.sh`: script to start the Docker environment.

# Local installation

Clone the repository and enter the directory:

```bash
git clone https://github.com/rprebet/albatross-2026
cd albatross-2026
```

## Without Docker

No root access is required.

### 1. Install Julia

Install Julia from the [Julia website](https://julialang.org/downloads/).

On Linux and macOS, Julia can also be installed with:

```bash
curl -fsSL https://install.julialang.org | sh
```

### 2. Install the required packages

<b>Do one and only one of the following.</b>

### Automatic installation

Run:

```bash
./install.sh
```


### Manual installation

Alternatively, start Julia:

```bash
julia
```

and run:

```julia
using Pkg
Pkg.add("AlgebraicSolving")
Pkg.add(["Plots", "IJulia"])

using AlgebraicSolving, Plots, IJulia
notebook()
```

The manual installation may cause version conflicts if packages are already installed in the global Julia environment.

### 3. Install Jupyter

If Julia asks whether to install Jupyter via Conda, press `Enter` (or answer `y`).

### 4. Start the practical session

Open `warmup.ipynb` and run the first cell. This performs some initial package precompilation.

If Jupyter is not working, the corresponding Julia files can be found in `julia_files/`.

# Docker

The Docker option provides a self-contained environment with Julia and the required packages.

Docker installation generally requires root access.

### 1. Install Docker

For example, on Ubuntu:

```bash
sudo apt install docker.io
```

or:

```bash
sudo snap install docker
```

### 2. Start the environment

Go to the directory where you want to download the practical-session files. You need write access to this directory.

Then run:

```bash
./install_docker.sh
```

By default, Jupyter is exposed on port `8888`. If this port is already in use, change `8888` to another available port in `install_docker.sh`.

### 3. Open Jupyter

The script will print a URL of the form:

```text
http://127.0.0.1:8888/tree?token=...
```

Open this URL in a browser. If necessary, copy the token from the URL and enter it manually.

### 4. Start the practical session

Open `warmup.ipynb` and run the first cell. This performs some initial package precompilation.

# GitHub Codespaces

The practical session can also be run online using GitHub Codespaces. A GitHub account is required.

1. Go to https://github.com/rprebet/albatross-2026.
2. Click **Code**, select **Codespaces**, then **Create codespace on main**.
3. Wait for the virtual machine and the required packages and extensions to be installed. The initial setup can take up to 10 minutes.
4. If prompted, click **Trust folder**.
5. Open `warmup.ipynb` and run the first cell. This performs some initial package precompilation.

After this step, the environment is ready for the practical session.
