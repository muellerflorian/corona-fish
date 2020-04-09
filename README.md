# corona-fish
Repository with ressources used for corona-fish project

## Organization

* `covfish`: python modules with helper functions
* `data`: contains different data-sets
  * `fasta`: fasta sequences against fish probes were designed
  * `genomes`: different genomes for local blast
* `workflows`: different analysis workflows.
    **TOOO**: workflows as jupyter notebooks.

## Getting started
The provided code uses **Python** and **R**.

* For Python, we recommend installing [**Miniconda** with Python](https://docs.conda.io/en/latest/miniconda.html):
choose Python 3.7 and your operating system.
* For R, we recommend installing [**R**](https://www.r-project.org/) and the free versio of [RStudio Desktop](https://rstudio.com/products/rstudio/download/). 

### Python 
Most of the provide code runs under python. 

We recommend creating a **dedicated environment** to run code in this analysis package.

To create an environment called `cov-fish`, open an anaconda prompt and type (Confirm with `y` when asked if you want to proceed (`Proceed ([y]/n)?`):
```
conda create --name cov-fish python=3.7
```

**Activate the environment**:
```
conda activate cov-fish
```

**Install code for this package**:

* From the web: `pip install git+git:https://github.com/muellerflorian/corona-fish`
* Or a local dev install: `pip install -e /path/to/package`


## Workflows

### Probe-design

The file `workflows\probe_design\probe-design-overview.html` contains a detailed description 
of the different steps in the analysis workflow.