# 📦 MASCDB - A global archive of MASC snowflake images and associated descriptors.

This repository provide the software to create, manipulate and analyze Multi Angle Snowflake Camera (MASC) data.

|                   |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Deployment        | [![PyPI](https://badge.fury.io/py/mascdb.svg?style=flat)](https://pypi.org/project/mascdb/) [![Conda](https://img.shields.io/conda/vn/conda-forge/mascdb.svg?logo=conda-forge&logoColor=white&style=flat)](https://anaconda.org/conda-forge/mascdb)                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Activity          | [![PyPI Downloads](https://img.shields.io/pypi/dm/mascdb.svg?label=PyPI%20downloads&style=flat)](https://pypi.org/project/mascdb/) [![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/pymascdb.svg?label=Conda%20downloads&style=flat)](https://anaconda.org/conda-forge/pymascdb)                                                                                                                                                                                                                                                                                                                                                                                   |
| Python Versions   | [![Python Versions](https://img.shields.io/badge/Python-3.10%20%203.11%20%203.12%20%203.13%20%203.14-blue?style=flat)](https://www.python.org/downloads/)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Supported Systems | [![Linux](https://img.shields.io/github/actions/workflow/status/ltelab/pymascdb/.github/workflows/tests.yml?label=Linux&style=flat)](https://github.com/ltelab/mascdb/actions/workflows/tests.yml) [![macOS](https://img.shields.io/github/actions/workflow/status/ltelab/pymascdb/.github/workflows/tests.yml?label=macOS&style=flat)](https://github.com/ltelab/pymascdb/actions/workflows/tests.yml) [![Windows](https://img.shields.io/github/actions/workflow/status/ltelab/mascdb/.github/workflows/tests_windows.yml?label=Windows&style=flat)](https://github.com/ltelab/mascdb/actions/workflows/tests_windows.yml)                                                      |
| Project Status    | [![Project Status](https://www.repostatus.org/badges/latest/active.svg?style=flat)](https://www.repostatus.org/#active)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| Build Status      | [![Tests](https://github.com/ltelab/pymascdb/actions/workflows/tests.yml/badge.svg?style=flat)](https://github.com/ltelab/pymascdb/actions/workflows/tests.yml) [![Lint](https://github.com/ltelab/dispymascdbdrodb/actions/workflows/lint.yml/badge.svg?style=flat)](https://github.com/ltelab/pymascdb/actions/workflows/lint.yml) [![Docs](https://readthedocs.org/projects/pymascdb/badge/?version=latest&style=flat)](https://pymascdb.readthedocs.io/en/latest/)                                                                                                                                                                                                            |
| Linting           | [![Black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat)](https://github.com/psf/black) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json&style=flat)](https://github.com/astral-sh/ruff) [![Codespell](https://img.shields.io/badge/Codespell-enabled-brightgreen?style=flat)](https://github.com/codespell-project/codespell)                                                                                                                                                                                                                                                     |
| Code Coverage     | [![Coveralls](https://coveralls.io/repos/github/ltelab/pymascdb/badge.svg?branch=main&style=flat)](https://coveralls.io/github/ltelab/pymascdb?branch=main) [![Codecov](https://codecov.io/gh/ltelab/pymascdb/branch/main/graph/badge.svg?style=flat)](https://codecov.io/gh/ltelab/pymascdb)                                                                                                                                                                                                                                                                                                                                                                                     |
| Code Quality      | [![Codefactor](https://www.codefactor.io/repository/github/ltelab/pymascdb/badge?style=flat)](https://www.codefactor.io/repository/github/ltelab/pymascdb) [![Codebeat](https://codebeat.co/badges/14ff831b-f064-4bdd-a2e2-72ffdf28a35a?style=flat)](https://codebeat.co/projects/github-com-ltelab-pymascdb-main) [![Codacy](https://app.codacy.com/project/badge/Grade/d823c50a7ad14268bd347b5aba384623?style=flat)](https://app.codacy.com/gh/ltelab/pymascdb/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade) [![Codescene](https://codescene.io/projects/36773/status-badges/code-health?style=flat)](https://codescene.io/projects/36773) |
| License           | [![License](https://img.shields.io/github/license/ltelab/pymascdb?style=flat)](https://github.com/ltelab/pymascdb/blob/main/LICENSE)              |

| Citation | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17886482.svg?style=flat)](https://doi.org10.5281/zenodo.7398284) |


[**Documentation**](https://pymascdb.readthedocs.io/en/latest/) | [**Data Archive**](https://doi.org/10.5281/zenodo.5578920)

MASC is an international joint effort collect and homogenize Multi Angle Snowflake Camera (MASC) data from around the world.

## 🛠️ Installation

### conda

pymascdb can be installed via [conda][conda_link] on Linux, Mac, and Windows.
Install the package by typing the following command in the terminal:

```bash
conda install -c conda-forge pymascdb
```

In case conda-forge is not set up for your system yet, see the easy to follow instructions on [conda-forge][conda_forge_link].

### pip

pymascdb can be installed also via [pip][pip_link] on Linux, Mac, and Windows.
On Windows you can install [WinPython][winpy_link] to get Python and pip running.

Then, install the pymascdb package by typing the following command in the terminal:

```bash
pip install pymascdb
```

To install the latest development version via pip, see the [documentation][dev_install_link].

### 📚 Tutorial

The folder `tutorials` provides code examples to explore the capabilities of pymascdb.
A selection of of jupyter notebooks illustrates a selection of pymascdb functionalities.
These jupyter notebooks tutorial can also be consulted in the [`online documentation`].

## 📖 Explore the MASCDB documentation

With this introduction, we just scratched the surface of the MASCDB software capabilities.
To discover more about the MASCDB products, the download, processing and analysis utilities, or how to contribute your own data to MASCDB,
please read the software documentation available at [https://pymascdb.readthedocs.io/en/latest/index.html](https://pymascdb.readthedocs.io/en/latest/index.html).

## 💭 Feedback and Contributing Guidelines

If you aim to contribute your data or discuss the future development of MASCDB,
feel free to also open a [GitHub Issue](https://github.com/ltelab/pymascdb/issues) or a
[GitHub Discussion](https://github.com/ltelab/pymascdb/discussions) specific to your questions or ideas.

## ✍️ Contributors

- [Gionata Ghiggi](https://people.epfl.ch/gionata.ghiggi)
- [Jacopo Grazioli](https://people.epfl.ch/jacopo.grazioli)
- [Alexis Berne](https://people.epfl.ch/alexis.berne?lang=en)

## Citation

You can cite the MASCDB project by:

> Gionata Ghiggi, Jacopo Grazioli, Alexis Berne (2025). ltelab/pymascdb Zenodo. https://doi.org10.5281/zenodo.7398284
> Grazioli, J., Ghiggi, G., & Berne, A. (2023). MASCDB, a database of images, descriptors and microphysical properties of individual snowflakes in free fall (1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.8046497
> Grazioli, J., Ghiggi, G., Billault-Roux, AC. et al. MASCDB, a database of images, descriptors and microphysical properties of individual snowflakes in free fall. Sci Data 9, 186 (2022). https://doi.org/10.1038/s41597-022-01269-7

If you want to cite a specific version of pymascdb, have a look at the [Zenodo software archive repository](https://zenodo.org/records/17886482).

## License

The content of this repository is released under the terms of the [MIT license](LICENSE).

[conda_forge_link]: https://github.com/conda-forge/mascdb-feedstock#installing-mascdb
[conda_link]: https://docs.conda.io/en/latest/miniconda.html
[dev_install_link]: https://pymascdb.readthedocs.io/en/latest/installation.html#installation-for-contributors
[pip_link]: https://pypi.org/project/mascdb
[winpy_link]: https://winpython.github.io/
[`online documentation`]: https://pymascdb.readthedocs.io/en/latest/examples.html
