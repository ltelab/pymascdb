# MASC_DB api - An API to query MASC data.

![Snowflake_SOM_Clustering](./figs/SOM_Clustering/MASC_SOM_Cluster.png)

The code in this repository provides an API to query, filter and visualize MASC data.
Documentation is available at https://pymascdb.readthedocs.io/en/latest/index.html

The data should be downloaded from Zenodo at: https://doi.org/10.5281/zenodo.5578920

ATTENTION: The code is subject to changes in the coming weeks / months.

The folder `tutorials` (will) provide line-by-line tutorials to explore the capabilities of the API of MASC_DB.
The tutorials are conceived for line-by-line execution for example using *spyder*. They cannot be run as scripts:

- Downloading MASCDB  [[`00_download_mascdb.py`]]
- Data manipulation and *api* functionalities [[`01_data_manipulation.py`]] 
- Exploratory data analysis (eda) [[`02_eda.py`]] 
- Image display and processing [[`03_image_processing.py`]] 

[`00_download_mascdb.py`]: https://github.com/ltelab/pymascdb/tree/master/tutorial/00_download_mascdb.py
[`01_data_manipulation.py`]: https://github.com/ltelab/pymascdb/tree/master/tutorial/01_data_manipulation.py
[`02_eda.py`]: https://github.com/ltelab/pymascdb/tree/master/tutorial/02_eda.py
[`03_image_processing.py`]: https://github.com/ltelab/pymascdb/tree/master/tutorial/03_image_processing.py


The folder `examples` provides a link to some jupyter notebooks with a minimal selection of functionalities, the same as the ones used in the [`online documentation`]

[`online documentation`]: https://pymascdb.readthedocs.io/en/latest/examples.html


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
