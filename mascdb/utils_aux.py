# -----------------------------------------------------------------------------.
# Copyright (c) 2021-2025 MASCDB developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the MIT License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License
# along with this program.  If not, see <https://opensource.org/license/mit/>.
# -----------------------------------------------------------------------------.
"""MASCDB auxiliary functions."""


####--------------------------------------------------------------------------.
#######################
#### Snowflake class ##
#######################
def get_snowflake_class_name_dict(method="Praz2017"):
    """
    Get snowflake class ID mapping from class name.

    Returns a dictionary mapping snowflake class names to their corresponding integer IDs
    according to the specified classification method.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017" based on
        https://amt.copernicus.org/articles/10/1335/2017/.

    Returns
    -------
    dict
        Dictionary mapping class names (str) to class IDs (int).
        For "Praz2017" method, includes:
        - "small_particle": 1
        - "columnar_crystal": 2
        - "planar_crystal": 3
        - "aggregate": 4
        - "graupel": 5
        - "columnar_planar_combination": 6

    Raises
    ------
    ValueError
        If the specified method is not available.

    """
    if method == "Praz2017":
        dict = {
            "small_particle": 1,
            "columnar_crystal": 2,
            "planar_crystal": 3,
            "aggregate": 4,
            "graupel": 5,
            "columnar_planar_combination": 6,
        }
    else:
        raise ValueError(f"Snowflake class dictionary not available for method {method}.")

    return dict


def get_snowflake_class_name_colors_dict(method="Praz2017"):
    """
    Get color mapping for snowflake class names.

    Returns a dictionary mapping snowflake class names to color names for visualization.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017".

    Returns
    -------
    dict
        Dictionary mapping class names (str) to color names (str).

    Raises
    ------
    ValueError
        If the specified method is not available.

    """
    if method == "Praz2017":
        dict = {
            "small_particle": "forestgreen",
            "columnar_crystal": "darkblue",
            "planar_crystal": "red",
            "aggregate": "orange",
            "graupel": "yellow",
            "columnar_planar_combination": "gray",
        }
    else:
        raise ValueError(f"Snowflake class dictionary not available for method {method}.")

    return dict


def get_snowflake_class_id_dict(method="Praz2017"):
    """
    Get snowflake class name mapping from class ID.

    Returns a dictionary mapping snowflake class IDs to their corresponding names.
    This is the inverse of get_snowflake_class_name_dict.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017".

    Returns
    -------
    dict
        Dictionary mapping class IDs (int) to class names (str).

    """
    dict = get_snowflake_class_name_dict(method=method)
    dict = {v: k for k, v in dict.items()}
    return dict


def get_snowflake_class_id_colors_dict(method="Praz2017"):
    """
    Get color mapping for snowflake class IDs.

    Returns a dictionary mapping snowflake class IDs to color names for visualization.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017".

    Returns
    -------
    dict
        Dictionary mapping class IDs (int) to color names (str).

    """
    colors_dict = get_snowflake_class_name_colors_dict(method=method)
    name_dict = get_snowflake_class_name_dict(method=method)
    dict = {name_dict[k]: v for k, v in colors_dict.items()}
    return dict


####--------------------------------------------------------------------------.
####################
#### Riming class ##
####################
def get_riming_class_name_dict(method="Praz2017"):
    """
    Get riming class ID mapping from class name.

    Returns a dictionary mapping riming class names to their corresponding integer IDs
    according to the specified classification method.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017" based on
        https://amt.copernicus.org/articles/10/1335/2017/.

    Returns
    -------
    dict
        Dictionary mapping class names (str) to class IDs (int).
        For "Praz2017" method, includes:
        - "undefined": 0
        - "unrimed": 1
        - "rimed": 2
        - "densely_rimed": 3
        - "graupel-like": 4
        - "graupel": 5

    Raises
    ------
    ValueError
        If the specified method is not available.

    """
    if method == "Praz2017":
        dict = {
            "undefined": 0,
            "unrimed": 1,
            "rimed": 2,
            "densely_rimed": 3,
            "graupel-like": 4,
            "graupel": 5,
        }
    else:
        raise ValueError(f"Riming class dictionary not available for method {method}.")

    return dict


def get_riming_class_name_colors_dict(method="Praz2017"):
    """
    Get color mapping for riming class names.

    Returns a dictionary mapping riming class names to color names for visualization.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017".

    Returns
    -------
    dict
        Dictionary mapping class names (str) to color names (str).

    Raises
    ------
    ValueError
        If the specified method is not available.

    """
    if method == "Praz2017":
        dict = {
            "undefined": "forestgreen",
            "unrimed": "darkblue",
            "rimed": "red",
            "densely_rimed": "orange",
            "graupel-like": "yellow",
            "graupel": "gray",
        }
    else:
        raise ValueError(f"Riming class dictionary not available for method {method}.")

    return dict


def get_riming_class_id_dict(method="Praz2017"):
    """
    Get riming class name mapping from class ID.

    Returns a dictionary mapping riming class IDs to their corresponding names.
    This is the inverse of get_riming_class_name_dict.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017".

    Returns
    -------
    dict
        Dictionary mapping class IDs (int) to class names (str).

    """
    dict = get_riming_class_name_dict(method=method)
    dict = {v: k for k, v in dict.items()}
    return dict


def get_riming_class_id_colors_dict(method="Praz2017"):
    """
    Get color mapping for riming class IDs.

    Returns a dictionary mapping riming class IDs to color names for visualization.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017".

    Returns
    -------
    dict
        Dictionary mapping class IDs (int) to color names (str).

    """
    colors_dict = get_riming_class_name_colors_dict(method=method)
    name_dict = get_riming_class_name_dict(method=method)
    dict = {name_dict[k]: v for k, v in colors_dict.items()}
    return dict


####--------------------------------------------------------------------------.
#####################
#### Melting class ##
#####################
def get_melting_class_name_dict(method="Praz2017"):
    """
    Get melting class ID mapping from class name.

    Returns a dictionary mapping melting class names to their corresponding integer IDs
    according to the specified classification method.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017".

    Returns
    -------
    dict
        Dictionary mapping class names (str) to class IDs (int).
        For "Praz2017" method, includes:
        - "dry": 0
        - "melting": 1

    Raises
    ------
    ValueError
        If the specified method is not available.

    """
    if method == "Praz2017":
        dict = {
            "dry": 0,
            "melting": 1,
        }
    else:
        raise ValueError(f"Melting class dictionary not available for method {method}.")

    return dict


def get_melting_class_name_colors_dict(method="Praz2017"):
    """
    Get color mapping for melting class names.

    Returns a dictionary mapping melting class names to color names for visualization.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017".

    Returns
    -------
    dict
        Dictionary mapping class names (str) to color names (str).

    Raises
    ------
    ValueError
        If the specified method is not available.

    """
    if method == "Praz2017":
        dict = {
            "dry": "darkblue",
            "melting": "orange",
        }
    else:
        raise ValueError(f"Melting class dictionary not available for method {method}.")

    return dict


def get_melting_class_id_dict(method="Praz2017"):
    """
    Get melting class name mapping from class ID.

    Returns a dictionary mapping melting class IDs to their corresponding names.
    This is the inverse of get_melting_class_name_dict.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017".

    Returns
    -------
    dict
        Dictionary mapping class IDs (int) to class names (str).

    """
    dict = get_melting_class_name_dict(method=method)
    dict = {v: k for k, v in dict.items()}
    return dict


def get_melting_class_id_colors_dict(method="Praz2017"):
    """
    Get color mapping for melting class IDs.

    Returns a dictionary mapping melting class IDs to color names for visualization.

    Parameters
    ----------
    method : str, optional
        Hydrometeor classification method. Default is "Praz2017".

    Returns
    -------
    dict
        Dictionary mapping class IDs (int) to color names (str).

    """
    colors_dict = get_melting_class_name_colors_dict(method=method)
    name_dict = get_melting_class_name_dict(method=method)
    dict = {name_dict[k]: v for k, v in colors_dict.items()}
    return dict


####--------------------------------------------------------------------------.
############################
#### Precipitation Class ###
############################
def get_precip_class_name_dict(method="Schaer2020"):
    """
    Get precipitation class ID mapping from class name.

    Returns a dictionary mapping precipitation class names to their corresponding integer IDs
    according to the specified classification method.

    Parameters
    ----------
    method : str, optional
        Precipitation classification method. Default is "Schaer2020".

    Returns
    -------
    dict
        Dictionary mapping class names (str) to class IDs (int).
        For "Schaer2020" method, includes:
        - "undefined": 0
        - "precip": 1
        - "mixed": 2
        - "blowing_snow": 3

    Raises
    ------
    ValueError
        If the specified method is not available.

    """
    if method == "Schaer2020":
        dict = {
            "undefined": 0,
            "precip": 1,
            "mixed": 2,
            "blowing_snow": 3,
        }
    else:
        raise ValueError(f"Precipitation class dictionary not available for method {method}.")

    return dict


def get_precip_class_name_colors_dict(method="Praz2017"):
    """
    Get color mapping for precipitation class names.

    Returns a dictionary mapping precipitation class names to color names for visualization.

    Parameters
    ----------
    method : str, optional
        Precipitation classification method. Default is "Praz2017".

    Returns
    -------
    dict
        Dictionary mapping class names (str) to color names (str).

    Raises
    ------
    ValueError
        If the specified method is not available.

    """
    if method == "Schaer2020":
        dict = {
            "undefined": "forestgreen",
            "precip": "darkblue",
            "mixed": "orange",
            "blowing_snow": "yellow",
        }
    else:
        raise ValueError(f"Precipitation class dictionary not available for method {method}.")

    return dict


def get_precip_class_id_dict(method="Schaer2020"):
    """
    Get precipitation class name mapping from class ID.

    Returns a dictionary mapping precipitation class IDs to their corresponding names.
    This is the inverse of get_precip_class_name_dict.

    Parameters
    ----------
    method : str, optional
        Precipitation classification method. Default is "Schaer2020".

    Returns
    -------
    dict
        Dictionary mapping class IDs (int) to class names (str).

    """
    dict = get_precip_class_name_dict(method=method)
    dict = {v: k for k, v in dict.items()}
    return dict


def get_precip_class_id_colors_dict(method="Schaer2020"):
    """
    Get color mapping for precipitation class IDs.

    Returns a dictionary mapping precipitation class IDs to color names for visualization.

    Parameters
    ----------
    method : str, optional
        Precipitation classification method. Default is "Schaer2020".

    Returns
    -------
    dict
        Dictionary mapping class IDs (int) to color names (str).

    """
    colors_dict = get_precip_class_name_colors_dict(method=method)
    name_dict = get_precip_class_name_dict(method=method)
    dict = {name_dict[k]: v for k, v in colors_dict.items()}
    return dict


####--------------------------------------------------------------------------.
#######################
#### Campaign Utils ###
#######################
def get_campaign_colors_dict():
    """Get a dictionary mapping campaign names to colors."""
    d_c = {
        "APRES3-2016": "forestgreen",
        "APRES3-2017": "lightgreen",
        "Davos-2015": "darkblue",
        "Davos-2019": "lightblue",
        "Valais-2016": "turquoise",
        "ICEPOP-2018": "violet",
        "PLATO-2019": "pink",
        "Jura-2019": "orange",
        "POPE-2020": "yellow",
        "ICEGENESIS-2021": "red",
    }
    return d_c


def get_campaign_names():
    """Get list of available campaign names."""
    return list(get_campaign_colors_dict().keys())


####--------------------------------------------------------------------------.
##############
#### Units ###
##############


def var_units():
    """
    Return dictionary containing units of MASC descriptors.

    Provides a comprehensive mapping of all MASC (Multi-Angle Snowflake Camera)
    descriptor variable names to their physical units.

    Returns
    -------
    dict
        Dictionary mapping variable names (str) to their units (str).
        Common units include: 'm' (meters), 'deg' (degrees), '-' (dimensionless),
        'pix' (pixels), 'class' (classification), 'boolean', etc.

    """
    units = {
        "datetime": "datetime",
        "index": "-",
        "flake_id": "-",
        "flake_number_tmp": "-",
        "pix_size": "m",
        "cam_id": "-",
        "quality_xhi": "-",
        "n_roi": "-",
        "flake_n_roi": "-",
        "area": "m**2",
        "perim": "m",
        "Dmean": "m",
        "Dmax": "m",
        "eq_radius": "m",
        "area_porous": "m**2",
        "area_porous_r": "-",
        "ell_fit_A": "m",
        "ell_fit_B": "m",
        "ell_fit_area": "m**2",
        "ell_fit_ori": "deg",
        "ell_fit_a_r": "-",
        "ell_fit_ecc": "-",
        "compactness": "-",
        "ell_in_A": "m",
        "ell_in_B": "m",
        "ell_in_area": "m**2",
        "ell_out_A": "m",
        "ell_out_B": "m",
        "ell_out_area": "m**2",
        "roundness": "-",
        "p_circ_out_r": "-",
        "rectangularity": "-",
        "bbox_width": "m",
        "bbox_len": "m",
        "rect_perim_ratio": "-",
        "rect_aspect_ratio": "-",
        "rect_eccentricity": "-",
        "solidity": "-",
        "convexity": "-",
        "hull_n_angles": "-",
        "p_circ_r": "-",
        "frac_dim_boxcounting": "-",
        "frac_dim_theoretical": "-",
        "nb_holes": "-",
        "skel_N_ends": "-",
        "skel_N_junc": "-",
        "skel_perim_ratio": "-",
        "skel_area_ratio": "pix**-1",
        "sym_P1": "-",
        "sym_P2": "-",
        "sym_P3": "-",
        "sym_P4": "-",
        "sym_P5": "-",
        "sym_P6": "-",
        "sym_Pmax_id": "-",
        "sym_P6_max_ratio": "-",
        "sym_mean": "pix",
        "sym_std": "pix",
        "sym_std_mean_ratio": "-",
        "intensity_mean": "-",
        "intensity_max": "-",
        "contrast": "-",
        "intensity_std": "-",
        "hist_entropy": "-",
        "local_std": "-",
        "local_intens": "-",
        "lap_energy": "-",
        "wavs": "-",
        "complexity": "-",
        "har_energy": "-",
        "har_contrast": "-",
        "har_corr": "-",
        "har_hom": "-",
        "roi_centroid_X": "pix",
        "roi_centroid_Y": "pix",
        "roi_width": "pix",
        "roi_height": "pix",
        "Dmax_ori": "deg",
        "Dmax_90": "m",
        "D90_r": "-",
        "riming_class_id": "class",
        "riming_class_prob": "-",
        "riming_class_name": "class string",
        "riming_deg_level": "-",
        "melting_class_id": "boolean",
        "melting_class_name": "class_string",
        "melting_prob": "-",
        "snowflake_class_name": "class string",
        "snowflake_class_id": "class",
        "snowflake_class_prob": "-",
        "flake_fallspeed": "m s**-1",
        "campaign": "-",
        "latitude": "deg_north",
        "longitude": "deg_east",
        "altitude": "m",
        "flake_quality_xhi": "-",
        "flake_Dmax": "m",
        "gan3d_mass": "kg",
        "gan3d_volume": "m**3",
        "gan3d_gyration": "m",
        "bs_normalized_angle": "-",
        "bs_mixing_ind": "-",
        "bs_precip_class_name": "class string",
        "bs_precip_class_id": "-",
        "env_T": "deg C",
        "env_P": "hPa",
        "env_DD": "deg",
        "env_FF": "m s**-1",
        "env_RH": "%",
        "hl_snowflake": "boolean",
        "hl_snowflake_class_id": "int",
        "hl_melting": "boolean",
        "hl_melting_class_id": "int",
        "hl_riming": "boolean",
        "hl_riming_class_id": "int",
    }
    return units


####--------------------------------------------------------------------------.
##############
#### Description ###
##############


def var_explanations():
    """
    Get dictionary containing verbose explanations of MASC descriptors.

    Provides detailed descriptions for all MASC (Multi-Angle Snowflake Camera)
    descriptor variables, including references to relevant publications.

    Returns
    -------
    dict
        Dictionary mapping variable names (str) to their detailed explanations (str).
        Explanations include physical meaning, calculation methods, and references
        to scientific publications where applicable.

    """
    explanations = {
        "datetime": "Datetime object of the measurement",
        "index": (
            "Index ranging from 0 to N, where N is the number of observations in the database. For unique identifications better is to use flake_id"  # noqa: E501
        ),
        "flake_id": (
            "Unique identifier of each measurement. It combines the datetime of measurement with the temporary internal flake number given by the MASC"  # noqa: E501
        ),
        "flake_number_tmp": "Temporary flake number. Incremental, but it resets upon reboot of the instrument. ",
        "pix_size": "Pixel size",
        "quality_xhi": (
            "Quality index of the ROI. Very good images above values of 9.  Reference is https://doi.org/10.5194/amt-10-1335-2017 (see Appendix B)"  # noqa: E501
        ),
        "cam_id": "ID of the CAM: 0, 1 or 2",
        "n_roi": (
            "Number of ROIs initially identified in the raw image of one camera. Note that all the processing downstream is referred to only one (the main) ROI"  # noqa: E501
        ),
        "flake_n_roi": "Average value of n_roi (see n_roi definition) over the three cameras ",
        "area": "ROI area. Descriptor 1 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)",
        "perim": "ROI perimeter. Descriptor 2 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)",
        "Dmean": (
            "ROI mean diameter. Mean value of x-width and y-height. Descriptor 3 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "Dmax": "ROI maximum dimension. Descriptor 4 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)",
        "eq_radius": (
            "ROI equi-areal radius. Radius of a circle having the same area of the ROI. Descriptor 5 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "area_porous": (
            "ROI area with internal holes removed. Descriptor 6 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "area_porous_r": (
            "Ratio between area_porous and area. Descriptor 7 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_fit_A": (
            "Major semidimension of the ellipse fitted on the ROI. Descriptor 8 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_fit_B": (
            "Minor semidimension of the ellipse fitted on the ROI. Descriptor 9 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_fit_area": (
            "Area of the ellipse fitted on the ROI. Descriptor 10 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_fit_ori": (
            "Orientation of the ellipse fitted on the ROI. Descriptor 11 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_fit_a_r": (
            "Axis ratio of the ellipse fitted on the ROI. Descriptor 12 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_fit_ecc": (
            "Eccentricity of the ellipse fitted on the ROI. Descriptor 13 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "compactness": (
            "Ratio between projected area and area of the fitted ellipse. Descriptor 14 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_in_A": (
            "Major semidimension of inscribed ellipse having same center and orientation as the fitted one. Descriptor 15 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_in_B": (
            "Minor semidimension of inscribed ellipse having same center and orientation as the fitted one. Descriptor 16 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_in_area": (
            "Area of inscribed ellipse having same center and orientation as the fitted one. Descriptor 17 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_out_A": (
            "Major semidimension of circumscribed ellipse having same center and orientation as the fitted one. Descriptor 18 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_out_B": (
            "Minor semidimension of circumscribed ellipse having same center and orientation as the fitted one. Descriptor 19 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "ell_out_area": (
            "Area of circumscribed ellipse having same center and orientation as the fitted one. Descriptor 20 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "roundness": (
            "Ratio between ROI area and area of circumscribed circle. Descriptor 30 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A) "  # noqa: E501
        ),
        "p_circ_out_r": (
            "Ratio between ROI perimeter and perimeter  of circumscribed circle. Descriptor 31 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A) "  # noqa: E501
        ),
        "rectangularity": (
            "Ratio between ROI area and area of bounding box. Descriptor 32 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "bbox_width": (
            "Width of bounding box of the ROI. Descriptor 33 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "bbox_len": (
            "Length of bounding box of the ROI. Descriptor 34 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "rect_perim_ratio": (
            "Ratio between ROI bounding box perimeter and ROI perimeter. Descriptor 35 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "rect_aspect_ratio": (
            "Aspect ratio of ROI bounding box. Descriptor 36 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "rect_eccentricity": (
            "Eccentricity of ROI bounding box. Descriptor 37 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "solidity": (
            "Ratio between ROI area and convex-hull area. Descriptor 38 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "convexity": (
            "Ratio between ROI perimetr and convex-hull perimeter. Descriptor 39 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "hull_n_angles": (
            "Number of vertices of the convex-hull of the ROI. Descriptor 40 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "p_circ_r": (
            "Ratio between ROI perimeter and perimeter of equivalent-area circle. Descriptor 41 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "frac_dim_boxcounting": (
            "Boxcounting ROI fractal dimension. Descriptor 42 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "frac_dim_theoretical": (
            "Theoretical fractal dimensions (calculated from area and perimeter). Descriptor 43 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "nb_holes": (
            "Number of holes identified inside the ROI (note that snowflake aggregates or crystals may have holes, it is not necessarily an artifact)."  # noqa: E501
        ),
        "skel_N_ends": (
            "Number of ends of the inner skeleton of the ROI. Descriptor 44 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "skel_N_junc": (
            "Number of junctions of the inner skeleton of the ROI. Descriptor 45 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "skel_perim_ratio": (
            "Ratio between skeleton length and ROI perimeter. Descriptor 46 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "skel_area_ratio": (
            "Ratio between skeleton length and ROI area. Descriptor 47 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "sym_P1": (
            "Standardized distance to centroid Fourier power spectrum comp. P1. Descriptor 49 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)  "  # noqa: E501
        ),
        "sym_P2": (
            "Standardized distance to centroid Fourier power spectrum comp. P2. Descriptor 50 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A) "  # noqa: E501
        ),
        "sym_P3": (
            "Standardized distance to centroid Fourier power spectrum comp. P3. Descriptor 51 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A) "  # noqa: E501
        ),
        "sym_P4": (
            "Standardized distance to centroid Fourier power spectrum comp. P4. Descriptor 52 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A) "  # noqa: E501
        ),
        "sym_P5": (
            "Standardized distance to centroid Fourier power spectrum comp. P5. Descriptor 53 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A) "  # noqa: E501
        ),
        "sym_P6": (
            "Standardized distance to centroid Fourier power spectrum comp. P6. Descriptor 54 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A) "  # noqa: E501
        ),
        "sym_Pmax_id": (
            "Maximum ID (1 to 6) of variables sym_P*. Descriptor 55 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "sym_P6_max_ratio": (
            "Ratio between sym_P6 and max(sym_P*). Descriptor 56 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "sym_mean": (
            "Mean distance to centroid. Descriptor 57 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"
        ),
        "sym_std": (
            "Standard deviation of distance to centroid. Descriptor 58 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "sym_std_mean_ratio": (
            "Ratio between sym_std and sym_mean. Descriptor 59 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "intensity_mean": (
            "ROI mean pixel brightness. Descriptor 60 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"
        ),
        "intensity_max": (
            "ROI maximum pixel brightness. Descriptor 61 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"
        ),
        "contrast": "ROI contrast. Descriptor 62 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A) ",
        "intensity_std": (
            "ROI standard deviation of pixel brightness. Descriptor 63 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "hist_entropy": (
            "Brightness histogram entropy. Descriptor 64 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"
        ),
        "local_std": (
            "Average greyscale ROI local standard deviation in a 3x3 window. Descriptor 65 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "local_intens": (
            "Average local intensity in a 3x3 window. Descriptor 66 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "lap_energy": (
            "Energy of the Laplacian. Descriptor 67 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A) "
        ),
        "wavs": (
            "Sum of wavelet coeffficent. Descriptor 68 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"
        ),
        "complexity": (
            "Particle complexity as in  https://doi.org/10.1002/2014GL061016. Descriptor 69 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)"  # noqa: E501
        ),
        "har_energy": "Haralick Energy. Descriptor 70 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)",
        "har_contrast": "Haralick contrast. Descriptor 71 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)",
        "har_corr": "Haralick correlation. Descriptor 72 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)",
        "har_hom": "Haralick homogeneity. Descriptor 73 of https://doi.org/10.5194/amt-10-1335-2017 (see Appendix A)",
        "roi_centroid_X": (
            "X-position of the centroid of the ROI within the ORIGINAL MASC image (before cropping around the ROI itself)"  # noqa: E501
        ),
        "roi_centroid_Y": (
            "Y-position of the centroid of the ROI within the ORIGINAL MASC image (before cropping around the ROI itself)"  # noqa: E501
        ),
        "roi_width": "X-size of the cropped ROI",
        "roi_height": "Y-size of the cropped ROI",
        "Dmax_ori": "Orientation of maximum dimension Dmax (Dmax includes points within the ROI)",
        "Dmax_90": (
            "Maximum dimension in the orthogonal direction of Dmax. See visual representation in https://doi.org/10.1175/JAM2398.1"
        ),
        "D90_r": "Axis ratio between Dmax and Dmax_90",
        "riming_class_id": (
            "ID of riming class (0 to 5) from https://doi.org/10.5194/amt-10-1335-2017. 0 meaning undefined, while the other classes can be found in the paper."  # noqa: E501
        ),
        "riming_class_prob": "Probability associated with riming_class_id / riming_class_name",
        "riming_class_name": (
            "Riming class name (undefined, unrimed, rimed, densely_rimed, graupel-like, graupel. As defined in https://doi.org/10.5194/amt-10-1335-2017"
        ),
        "riming_deg_level": (
            "Continuously varying riming degree level (from 0, unrimed to 1 graupel). Variable named R_c in https://doi.org/10.5194/amt-10-1335-2017"
        ),
        "melting_class_id": (
            "ID of melting class (0: not melting, 1: melting). From the method of https://doi.org/10.5194/amt-10-1335-2017"
        ),
        "melting_class_name": (
            "Melting class name: (dry vs melting). From the method of https://doi.org/10.5194/amt-10-1335-2017"
        ),
        "melting_prob": (
            "Probability that a given particle or triplet is melting according to the method of https://doi.org/10.5194/amt-10-1335-2017. If rounded, this variable corresponds to melting_class_id   "  # noqa: E501
        ),
        "snowflake_class_name": (
            "Name of hydrometeor class associated to the ROI or to the triplet of ROIs according to the method of https://doi.org/10.5194/amt-10-1335-2017. (small_particle, columnar_crystal, planar_crystal, aggregate, graupel, columnar_planar_combination)"  # noqa: E501
        ),
        "snowflake_class_id": (
            "ID of hydrometeor class associated to the ROI or to the triplet of ROIs according to the method of https://doi.org/10.5194/amt-10-1335-2017. (small_particle: 1, columnar_crystal: 2, planar_crystal: 3, aggregate: 4, graupel: 5, columnar_planar_combination: 6)"  # noqa: E501
        ),
        "snowflake_class_prob": "Probability associated with snowflake_class_id  / snowflake_class_name",
        "flake_fallspeed": "Fall speed as recorded by the MASC infrared sensors.",
        "campaign": "String indicating the name of the masurement campaign where the MASC was deployed",
        "latitude": "WGS84 latitude",
        "longitude": "WGS84 longitude",
        "altitude": "Altitude on mean sea level",
        "flake_quality_xhi": (
            "Quality index of a triplet (mean value over the 3 ROIs, one for each camera). Very good images above values of 9.  Reference is https://doi.org/10.5194/amt-10-1335-2017 (see Appendix B)"  # noqa: E501
        ),
        "flake_Dmax": (
            "Maximum value of Dmax among the three individual ROIs (one for each camera). Used as a proxy for the true Dmax of the snowflake"  # noqa: E501
        ),
        "gan3d_mass": (
            "Estimated mass of the snowflake (when an estimation is possible). It is an output from the method of https://doi.org/10.5194/amt-2021-176 "  # noqa: E501
        ),
        "gan3d_volume": (
            "Estimated 3D volume of the snowflake (when an estimation is possible). It is an output from the method of https://doi.org/10.5194/amt-2021-176"
        ),
        "gan3d_gyration": (
            "Estimated radius of gyration of the snowflake (when an estimation is possible). It is an output from the method of https://doi.org/10.5194/amt-2021-176"  # noqa: E501
        ),
        "bs_normalized_angle": (
            "Blowing snow normalized angle. It is a parameter described in https://doi.org/10.5194/tc-14-367-2020 to discriminate between snow and blowing snow environments. If it is lower than 0.193 it indicates precipitation. Above 0.881 blowing snow. Mixed environments in between those values."  # noqa: E501
        ),
        "bs_mixing_ind": (
            "Blowing snow mixing index. It is a parameter described in https://doi.org/10.5194/tc-14-367-2020 . Defined only in case the method predicts a mix of precipitation and blowing snow. It ranges from 0 (precipitation) to 1 (pure blowing snow)"  # noqa: E501
        ),
        "bs_precip_class_name": (
            "Blowing snow precipitation class (undefined, precip, mixed, blowing_snow). Reference: https://doi.org/10.5194/tc-14-367-2020 "  # noqa: E501
        ),
        "bs_precip_class_id": (
            "Blowing snow precipitation ID (0: undefined, 1: precip, 2: mixed, 3: blowing_snow). Reference: https://doi.org/10.5194/tc-14-367-2020  "  # noqa: E501
        ),
        "env_T": "Environmental temperature in the proximity of the instrument",
        "env_P": "Environmental pressure in the proximity of the instrument",
        "env_DD": "Wind direction in the proximity of the instruments",
        "env_FF": "Wind speed (minute or minutes scale) in the proximity of the instrument  ",
        "env_RH": "Relative Humidity in the proximity of the instrument",
        "hl_snowflake": (
            "Boolean flag indicating if the ROI of a given cam was part of the hydrometeor classification human label (HL) trainingset of https://doi.org/10.5194/amt-10-1335-2017 "  # noqa: E501
        ),
        "hl_snowflake_class_id": (
            "Human label (HL) snowflake_class_id used in the training set of https://doi.org/10.5194/amt-10-1335-2017 "
        ),
        "hl_melting": (
            "Boolean flag indicating if the ROI of a given cam was part of the melting identification human label (HL) training set of https://doi.org/10.5194/amt-10-1335-2017"  # noqa: E501
        ),
        "hl_riming": (
            "Boolean flag indicating if the ROI of a given cam was part of the riming classification human label (HL) training set of https://doi.org/10.5194/amt-10-1335-2017"  # noqa: E501
        ),
        "hl_riming_class_id": (
            "Human label (HL) riming_class_id used in the training set of https://doi.org/10.5194/amt-10-1335-2017 "
        ),
    }
    return explanations


####--------------------------------------------------------------------------.
###########################
#### Dataframe columns  ###
###########################


def get_vars_gan3d():
    """
    Retrieve the list of all GAN3D variables.

    Returns variable names related to 3D mass, volume, and gyration estimates
    from the GAN3D method.

    Returns
    -------
    list of str
        List of GAN3D variable names: ['gan3d_mass', 'gan3d_volume', 'gan3d_gyration'].

    References
    ----------
    https://doi.org/10.5194/amt-2021-176

    """
    variables = [
        "gan3d_mass",
        "gan3d_volume",
        "gan3d_gyration",
    ]
    return variables


def get_vars_env():
    """
    Retrieve the list of all environmental variables.

    Returns variable names for meteorological measurements in the proximity
    of the MASC instrument.

    Returns
    -------
    list of str
        List of environmental variable names: temperature, pressure, wind direction,
        wind speed, and relative humidity.

    """
    variables = [
        "env_T",
        "env_P",
        "env_DD",
        "env_FF",
        "env_RH",
    ]
    return variables


def get_vars_blowing_snow():
    """
    Retrieve the list of all blowing snow variables.

    Returns variable names related to blowing snow classification and mixing indices.

    Returns
    -------
    list of str
        List of blowing snow variable names including normalized angle, mixing index,
        and precipitation class identifiers.

    References
    ----------
    https://doi.org/10.5194/tc-14-367-2020

    """
    variables = [
        "bs_normalized_angle",
        "bs_mixing_ind",
        "bs_precip_class_name",
        "bs_precip_class_id",
    ]
    return variables


def get_vars_location():
    """
    Retrieve the list of all location variables.

    Returns variable names for spatiotemporal information about measurements.

    Returns
    -------
    list of str
        List of location variable names: datetime, campaign, latitude, longitude,
        and altitude.

    """
    variables = [
        "datetime",
        "campaign",
        "latitude",
        "longitude",
        "altitude",
    ]
    return variables


def get_vars_class():
    """
    Retrieve the list of all classification variables.

    Returns variable names for snowflake classification including riming, melting,
    and hydrometeor type classifications with associated probabilities.

    Returns
    -------
    list of str
        List of classification variable names including class IDs, names,
        probabilities, and riming degree levels.

    """
    variables = [
        "riming_class_name",
        "riming_class_id",
        "riming_class_prob",
        "riming_deg_level",
        "melting_class_id",
        "melting_class_name",
        "melting_prob",
        "snowflake_class_name",
        "snowflake_class_id",
        "snowflake_class_prob",
    ]
    return variables


def get_vars_class_ids():
    """
    Retrieve the list of class ID variables.

    Returns variable names for integer class identifiers.

    Returns
    -------
    list of str
        List of class ID variable names: snowflake_class_id, riming_class_id,
        and melting_class_id.

    """
    variables = [
        "snowflake_class_id",
        "riming_class_id",
        "melting_class_id",
    ]
    return variables


def get_vars_class_names():
    """
    Retrieve the list of class name variables.

    Returns variable names for string class labels.

    Returns
    -------
    list of str
        List of class name variable names: snowflake_class_name, riming_class_name,
        and melting_class_name.

    """
    variables = [
        "snowflake_class_name",
        "riming_class_name",
        "melting_class_name",
    ]
    return variables


def get_vars_cam_info():
    """
    Retrieve the list of MASC camera information variables.

    Returns variable names for metadata about individual camera captures.

    Returns
    -------
    list of str
        List of camera information variable names including index, datetime,
        flake ID, pixel size, camera ID, and quality index.

    """
    variables = [
        "index",
        "datetime",
        "flake_id",
        "flake_number_tmp",
        "pix_size",
        "cam_id",
        "quality_xhi",
        # event_id
        # event_duration
    ]
    return variables


def get_vars_cam_descriptors():
    """
    Retrieve the list of MASC camera descriptors.

    Returns variable names for geometric, morphological, and textural descriptors
    computed from individual camera ROI (Region of Interest) images.

    Returns
    -------
    list of str
        Comprehensive list of descriptor variable names including shape features,
        texture features, symmetry features, and complexity measures.

    References
    ----------
    Descriptors are detailed in Appendix A of:
    https://doi.org/10.5194/amt-10-1335-2017

    """
    # TODO: Or to infer from self.cam0.columns - get_vars_cam_info() - get_vars_class()
    # --> So that works if people add stuff ...
    variables = [
        "n_roi",
        "area",
        "perim",
        "Dmean",
        "Dmax",
        "eq_radius",
        "area_porous",
        "area_porous_r",
        "ell_fit_A",
        "ell_fit_B",
        "ell_fit_area",
        "ell_fit_ori",
        "ell_fit_a_r",
        "ell_fit_ecc",
        "compactness",
        "ell_in_A",
        "ell_in_B",
        "ell_in_area",
        "ell_out_A",
        "ell_out_B",
        "ell_out_area",
        "roundness",
        "p_circ_out_r",
        "rectangularity",
        "bbox_width",
        "bbox_len",
        "rect_perim_ratio",
        "rect_aspect_ratio",
        "rect_eccentricity",
        "solidity",
        "convexity",
        "hull_n_angles",
        "p_circ_r",
        "frac_dim_boxcounting",
        "frac_dim_theoretical",
        "nb_holes",
        "skel_N_ends",
        "skel_N_junc",
        "skel_perim_ratio",
        "skel_area_ratio",
        "sym_P1",
        "sym_P2",
        "sym_P3",
        "sym_P4",
        "sym_P5",
        "sym_P6",
        "sym_Pmax_id",
        "sym_P6_max_ratio",
        "sym_mean",
        "sym_std",
        "sym_std_mean_ratio",
        "intensity_mean",
        "intensity_max",
        "contrast",
        "intensity_std",
        "hist_entropy",
        "local_std",
        "local_intens",
        "lap_energy",
        "wavs",
        "complexity",
        "har_energy",
        "har_contrast",
        "har_corr",
        "har_hom",
        "roi_centroid_X",
        "roi_centroid_Y",
        "roi_width",
        "roi_height",
        "Dmax_ori",
        "Dmax_90",
        "D90_r",
    ]
    return variables
