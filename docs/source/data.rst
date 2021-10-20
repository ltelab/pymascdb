.. _data:

Data 
=======================================

*mascdb* is built as an API to manipulate and operate with 
the dataset available at the following link (TODO: Zenodo link)

The dataset includes four *parquet* files, where scalar descriptors
of snowflake images are stored and a *Zarr* (zipped, to unzip) 
storage where the actual grayscale image triplets of snowflakes
in free fall are stored. 


.. toctree::
   :maxdepth: 2
   
   triplet
   cam
   zarr
