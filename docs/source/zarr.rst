.. _zarr:

MASCDB Zarr store
=======================================
The *Zarr* store *MASCdb.zarr* contains the grayscale images (0: black to 255 white) of triplets
of snowflake views, as recorded by the Multi-Angle-Snowflake-Camera.
The dimension of the *Zarr* store are *N_obs x 1024 x 1024 x 3*: for each triplet  in the database
and for each of the 3 views, an image of fixed size *1024 x 1024* is given, with the identified ROI
roughly centered in the middle of each image.
