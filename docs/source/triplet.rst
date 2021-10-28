.. _triplet:

MASCDB triplet data
=======================================
The file *MASCdb_triplet.parquet* contains the attributes listed in the table below.

+----------------------+-------+----------+-----------------------+-------------------------------------------------+
|    Parameter         | Units |     Type |     Long name         |     Reference / Format / Algorithm              |
+======================+=======+==========+=======================+=================================================+
|           **Global information**                                                                                  |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| datetime             |       | datetime |                       |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| campaign             |       | string   | Field campaign string |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| latitude             | deg   | float    | WGS84 latitude        |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| longitude            | deg   | float    | WGS84 longitude       |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| altitude             | m     | float    |                       |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
|            **Flake information**                                                                                  |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| flake_id             |       | string   | | Unique flake ID     | | e.g. 2015.02.10_11.55.10_flake_4              |
|                      |       |          | |                     | | YYYY.MM.DD_HH.mm.ss_flake_flake_number_tmp    |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| flake_number_tmp     |       | string   | | Temporary flake ID  |                                                 |
|                      |       |          | | (not unique)        |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| flake_quality_xhi    |       | float    | | Average quality ind | "xi" in Praz et al, 2017                        |
|                      |       |          | | (on the three cams) |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| flake_fallspeed      | m/s   | float    | Recorded fallspeed    |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| flake_n_roi          |       | int      | | Average N. of ROIs  |                                                 |
|                      |       |          | | (on the three cams) |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| flake_Dmax           | m     | float    | | Maximum Dmax        | Table A1:4 Praz et al, 2017                     |
|                      |       |          | | (on the three cams) |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
|      **Riming estimation information**                                                                            |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| riming_deg_level     |       | float    | | Continuous riming   | R_c in Praz et al, 2017                         |
|                      |       |          | | degree level        |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| riming_class_id      |       | int      | | Discrete riming     | | Praz et al, 2017                              |
|                      |       |          | | degree class ID     | | 0: undefined, 1: unrimed, 2: rimed            |
|                      |       |          | |                     | | 3: densely-rimed, 4: graupel-like, 5: graupel |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| riming_class_prob    |       | float    | | Riming classif.     | Praz et al, 2017                                |
|                      |       |          | | probability         |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| riming_class_name    |       | string   | | Discrete riming     | See riming_class_id                             |
|                      |       |          | | degree class name   |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
|      **Melting estimation information**                                                                           |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| melting_class_id     |       | int      | | Discrete melting    | | Praz et al, 2017                              |
|                      |       |          | | class ID            | | 0: dry, 1: melting                            |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| melting_prob         |       | float    |   Melting probability | | Praz et al, 2017                              |
|                      |       |          |                       | | If rounded, it yields melting_class_id        |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| melting_class_name   |       | string   | | Discrete melting    | See melting_class_id                            |
|                      |       |          | | class name          |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| **Hydrometeor type estimation information**                                                                       |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| snowflake_class_name |       | string   | | Hydrometeor         | | Praz et al, 2017 citePraz                     |
|                      |       |          | | class name          | | See snowflake_class_id                        |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| snowflake_class_id   |       | int      | | Hydrometeor         | | Praz et al, 2017 citePraz                     |
|                      |       |          | | class ID            | | 1: small_particle, 2: columnar_crystal,       |
|                      |       |          |                       | | 3: planar_crystal, 4: aggregate,              |
|                      |       |          |                       | | 5: graupel, 6: columnar_planar_combination    |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| snowflake_class_prob |       | float    | | Classification      |                                                 |
|                      |       |          | | probability         |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
|   **3D reconstruction / mass estimation**                                                                         |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| gan3d_mass           | kg    | float    | Estimated mass        | Leinonen et al, 2021                            |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| gan3d_volume         | m^3   | float    | Estimated volume      | Leinonen et al, 2021                            |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| gan3d_gyration       | m     | float    | Estimated gyration rd | Leinonen et al, 2021                            |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
|   **Co-located environmental information**                                                                        |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| env_T                | deg C | float    | Air temperature       |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| env_P                | hPa   | float    | Pressure              |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| env_DD               | deg   | float    | | Wind direction      |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| env_FF               | m/s   | float    | Wind speed            |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| env_RH               | \%    | float    | Relative humidity     |                                                 |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
|         **Blowing snow estimation**                                                                               |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| bs_normalized_angle  |       | float    | | Blowing Snow        | | Schaer et al 2020                             |
|                      |       |          | | normalized angle    | | Pure precip. if < 0.193, Pure BS if > 0.881   |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| bs_mixing_ind        |       | float    | | Blowing snow        | | Schaer et al 2020                             |
|                      |       |          | | mixing index        | | Only defined in mixed BS/precip environments  |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| bs_precip_class_name |       | string   | | Blowing snow        | | Schaer et al 2020                             |
|                      |       |          | | class name          | | See bs_precip_class_id                        |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+
| bs_precip_class_id   |       | int      | | Blowing snow        | | Schaer et al 2020,                            |
|                      |       |          | | class ID            | | 0: undefined, 1: precip, 2: mixed,            |
|                      |       |          |                       | | 3: blowing_snow                               |
+----------------------+-------+----------+-----------------------+-------------------------------------------------+

References
=========================================

- `Praz et al, 2017 <https://zenodo.org/record/5578921#.YXqUeJuxVH4>`_: Praz, C., Roulet, Y.-A., and Berne, A.: Solid hydrometeor classification and riming degree estimation from pictures collected with a Multi-Angle Snowflake Camera, Atmos. Meas. Tech., 10, 1335–1357, https://doi.org/10.5194/amt-10-1335-2017, 2017.

- `Schaer et al, 2020 <https://tc.copernicus.org/articles/14/367/2020/>`_: Schaer, M., Praz, C., and Berne, A.: Identification of blowing snow particles in images from a Multi-Angle Snowflake Camera, The Cryosphere, 14, 367–384, https://doi.org/10.5194/tc-14-367-2020, 2020. 

- `Leinonen et al, 2021 <https://amt.copernicus.org/articles/14/6851/2021/amt-14-6851-2021.html>`_: Leinonen, J., Grazioli, J., and Berne, A.: Reconstruction of the mass and geometry of snowfall particles from multi-angle snowflake camera (MASC) images, Atmos. Meas. Tech., 14, 6851–6866, https://doi.org/10.5194/amt-14-6851-2021, 2021. 


