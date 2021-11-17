# MASC-DB API - An API to query MASC data.

![Snowflake_SOM_Clustering](./figs/SOM_Clustering/MASC_SOM_Cluster.png)

The code in this repository provides an API to query, filter and visualize MASC data.
Documentation is available at https://pymascdb.readthedocs.io/en/latest/index.html

ATTENTION: The code is subject to changes in the coming weeks / months.

The folder `tutorials` (will) provide jupyter notebooks describing various features of MASC-DB api.
- Downloading MASC DB  [[`download.ipynb`]]
- Read and filter MASC DB [[`read_and_filter.ipynb`]]
- Exploratory Data Analysis of MASC DB [[`eda.ipynb`]]

[`download.ipynb`]: https://nbviewer.jupyter.org/github/deepsphere/deepsphere-weather/blob/outputs/tutorials/spherical_grids.ipynb
[`eda.ipynb`]: https://nbviewer.jupyter.org/github/deepsphere/deepsphere-weather/blob/outputs/tutorials/interpolation_pooling.ipynb

The folder `experiments` (will) provide examples for:
- Latent Space Exploration with UMAP / PCA
- Snowflake clustering with Self-Organizing Maps
- Snowflake representation learning using Info-GAN 
- 3D snowflake mass reconstruction 
- Multinomial logistic classification
- CNN classification

## Installation

For a local installation, follow the below instructions.

1. Clone this repository.
   ```sh
   git clone https://github.com/ltelab/pymascdb.git
   cd pymascdb
   ```

2. Install the dependencies using conda:
   ```sh
   conda env create -f environment.yml
   ```
3. Activate the mascdb conda environment 
   ```sh
   conda activate mascdb
   ```
   
4. Just for info... to update the environment.yml: 
   ```sh
   conda env export > environment.yml
   ```

5. With *mascdb* environment activated, install the package:
   ```sh
   python setup.py install
   ```

In this way, *mascdb* can simply be loaded, from any directory, as:
   ```python
   import mascdb.api
   from mascdb.api import MASC_DB
   ```

## References 

- [Manuscript](https://XXXX)
- [Slides](https://XXXX)
- [Presentation](https://XXXX)

## Contributors

* [Jacopo Grazioli](https://people.epfl.ch/jacopo.grazioli) 
* [Gionata Ghiggi](https://people.epfl.ch/gionata.ghiggi)

## License

The content of this repository is released under the terms of the [MIT license](LICENSE.txt).
