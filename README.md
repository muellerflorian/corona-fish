# corona-fish

This repository provide the code used to design smiFISH probes against SARS-cov-2., described in our paper, as well as up to date protocols:

__Title:__ to be defined.

__Authors:__ to be established.

## Organization of repository

* `covfish`: python modules with helper functions.
* `data`: contains data-sets required for the probe-design.
  * `fasta`: fasta sequences against fish probes were designed
  * `genomes`: different genomes for local blast
* `docs`: folder containing the documentation.
* `workflows`: different analysis workflows.

## Getting started

The provided code uses **Python** and **R**.

* For Python, we recommend installing [**Miniconda** with Python](https://docs.conda.io/en/latest/miniconda.html):
choose Python 3.7 and your operating system.
* For R, we recommend installing [**R**](https://www.r-project.org/) and the free versio of [RStudio Desktop](https://rstudio.com/products/rstudio/download/).

### Creating a dedicated Python environment

Most of the provide code runs under Python.

We recommend creating a **dedicated environment** to run code in this analysis package.

To create an environment called `cov-fish`, open an anaconda prompt and type (Confirm with `y` when asked if you want to proceed (`Proceed ([y]/n)?`):

``` Shell
conda create --name cov-fish python=3.7
```

**Activate the environment**:

``` Shell
conda activate cov-fish
```

**Install code for this package**:

* Pip install from the web:
  
  ``` Shell
  pip install git+git:https://github.com/muellerflorian/corona-fish
  ```

* Or a local pip install for development:
  
  ``` Shell
  pip install -e /path/to/package`
  ```

## Workflow for probe design

For a detailed description of the design process for the probes, please consult the dedicated 
[documentation](docs/probe-design-overview.md).

## Protocols

The folder `protocols` contains tried and tested protocols for smiFISH with the SARS-Cov-2 probes.