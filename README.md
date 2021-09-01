# MASC-DB API - An API to query MASC data.

![Snowflake_Latent_Interpolation](./figs/Snowflake_Latent_Interpolation.gif)

The code in this repository provides an API to query, filter and visualize MASC data.

ATTENTION: The code is subject to changes in the coming weeks / months.

The folder `tutorials` (will) provide jupyter notebooks describing various features of MASC-DB api.
- Downloading MASC DB  [[`download.ipynb`]]
- Read and filter MASC DB [[`read_and_filter.ipynb`]]
- Exploratory Data Analysis of MASC DB [[`eda.ipynb`]]

[`download.ipynb`]: https://nbviewer.jupyter.org/github/deepsphere/deepsphere-weather/blob/outputs/tutorials/spherical_grids.ipynb
[`eda.ipynb`]: https://nbviewer.jupyter.org/github/deepsphere/deepsphere-weather/blob/outputs/tutorials/interpolation_pooling.ipynb

The folder `experiments` (will) provide examples for:
- Multinomial logistic classification
- CNN classification
- Info-GAN clustering
- Latent Space Exploration with UMAP / PCA
- 3D snowflake mass reconstruction 

## Installation

For a local installation, follow the below instructions.

1. Clone this repository.
   ```sh
   git clone https://github.com/xxxx/xxxx.git
   cd xxxxx
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

## References 

- [Manuscript](https://XXXX)
- [Slides](https://XXXX)
- [Presentation](https://XXXX)

## Contributors

* [Jacopo Grazioli](https://people.epfl.ch/jacopo.grazioli) 
* [Gionata Ghiggi](https://people.epfl.ch/gionata.ghiggi)

## License

The content of this repository is released under the terms of the [MIT license](LICENSE.txt).
