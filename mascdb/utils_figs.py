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
"""MASCDB Visualization Utilities."""
# -----------------------------------------------------------------------------.
import matplotlib as mpl
import numpy as np
import pandas as pd


def cm2inch(*tupl):
    """Convert centimeters to inches.

    Parameters
    ----------
    *tupl : float or tuple of float
        Dimensions in centimeters. Can be individual numbers or a single tuple.

    Returns
    -------
    tuple of float
        Dimensions converted to inches.

    Examples
    --------
    >>> cm2inch(10, 20)
    (3.937, 7.874)
    >>> cm2inch((10, 20))
    (3.937, 7.874)
    """
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    return tuple(i / inch for i in tupl)


def minmax(x):
    """Return the minimum and maximum values of an array.

    Parameters
    ----------
    x : array-like
        Input array.

    Returns
    -------
    list of float
        List containing [minimum, maximum] values.
    """
    return [np.min(x), np.max(x)]


def get_c_cmap_from_color_dict(color_dict, labels):
    """Create color mapping and colormap for scatter plots from a color dictionary.

    This function converts a dictionary of label-to-color mappings into integer color
    indices and a matplotlib colormap suitable for use with plt.scatter.

    Parameters
    ----------
    color_dict : dict
        Dictionary mapping labels to color names or hex values.
    labels : array-like
        Array of labels corresponding to data points.

    Returns
    -------
    list
        A list containing [c, cmap] where:

            - c : numpy.ndarray
                Integer array of color indices for each label.
            - cmap : matplotlib.colors.ListedColormap
                Colormap object with unique colors from the dictionary.

    Examples
    --------
    >>> color_dict = {"A": "red", "B": "blue", "C": "green"}
    >>> labels = ["A", "B", "A", "C"]
    >>> c, cmap = get_c_cmap_from_color_dict(color_dict, labels)

    """
    c_names = [color_dict[x] for x in labels]
    # Retrieve c integer values
    c, c_unique_name = pd.factorize(c_names, sort=False)
    # Create cmap
    cmap = mpl.colors.ListedColormap(c_unique_name)
    # Return object
    return [c, cmap]


def get_legend_handles_from_colors_dict(colors_dict, marker="o"):
    """Create legend handles from a color dictionary for matplotlib legends.

    This function generates legend handles that can be used with matplotlib.pyplot.legend
    to create custom legend entries based on a color dictionary.

    Parameters
    ----------
    colors_dict : dict
        Dictionary mapping labels to colors (color names or hex values).
    marker : str, optional
        Marker style for legend entries. Default is "o" (filled circle).
        Options include:
        - "o" : filled circle
        - "s" : filled square
        - "PATCH" : filled large rectangle
        - Any valid matplotlib marker

    Returns
    -------
    list
        List of matplotlib handle objects (Line2D or Patch) for use in legend.

    Examples
    --------
    >>> colors = {"Category A": "red", "Category B": "blue"}
    >>> handles = get_legend_handles_from_colors_dict(colors, marker="s")
    >>> plt.legend(handles=handles)
    """
    import matplotlib as mpl

    if marker == "PATCH":
        # PATCH ('filled large rectangle')
        handles = []
        for key in colors_dict:
            data_key = mpl.patches.Patch(facecolor=colors_dict[key], edgecolor=colors_dict[key], label=key)
            handles.append(data_key)
    else:
        # Classical Markers
        handles = []
        for key in colors_dict:
            data_key = mpl.lines.Line2D(
                [0],
                [0],
                linewidth=0,
                marker=marker,
                label=key,
                markerfacecolor=colors_dict[key],
                markeredgecolor=colors_dict[key],
                markersize=3,
            )
            handles.append(data_key)
    return handles


def get_colors_from_cmap(x, cmap_name="Spectral", vmin=None, vmax=None, nan_color=None):
    """Map numeric values to colors using a matplotlib colormap.

    This function converts numeric values to hexadecimal color codes based on a specified
    colormap. It handles both arrays and dictionaries, and provides special handling for
    NaN values.

    Parameters
    ----------
    x : array-like, dict, or float
        Numeric values to map to colors. Can be a single value, array, or dictionary
        where values are numeric.
    cmap_name : str, optional
        Name of the matplotlib colormap to use. Default is "Spectral".
        Examples: "Spectral", "viridis", "plasma", "Wistia", "Green".
    vmin : float, optional
        Minimum value for colormap normalization. Default is None.
    vmax : float, optional
        Maximum value for colormap normalization. Default is None.
    nan_color : str, optional
        Hexadecimal color code to use for NaN values. If None, NaN values are
        marked as "NaN" string. Default is None.

    Returns
    -------
    numpy.ndarray or dict
        Hexadecimal color codes corresponding to input values. Returns a dictionary
        if input was a dictionary, otherwise returns an array.

    Examples
    --------
    >>> get_colors_from_cmap(2, cmap_name="Spectral", vmin=0, vmax=8)
    '#...'  # hex color
    >>> get_colors_from_cmap([0, 2, 3, 5], cmap_name="Wistia", vmin=0, vmax=5)
    array(['#...', '#...', '#...', '#...'], dtype='<U7')
    >>> get_colors_from_cmap({"a": 1, "b": 5}, cmap_name="viridis", vmin=0, vmax=10)
    {'a': '#...', 'b': '#...'}

    Notes
    -----
    NaN values in the input are preserved and can be assigned a custom color using
    the nan_color parameter.
    """
    # # - mpl.cm.<Spectral> : colormaps are defined by default between 0 and 255  mpl.cm.Spectral(257)
    flag_dict = False
    # Preprocess x if dictionary
    if isinstance(x, dict):
        flag_dict = True
        keys = list(x.keys())
        x = np.asarray(list(x.values()))
    # Get index with NaN
    idx_nan = np.isnan(x)
    # Retrieve colormap
    cmap = mpl.cm.get_cmap(cmap_name)
    # Rescale x and assign colormap values
    rgb_val = cmap((x - vmin) / vmax)
    # norm_fun = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    # rgb_val = cmap(norm_fun(x))
    # ----------------------------------------------------------------.
    ######################
    ## Convert to hex ####
    ######################
    # If only a single value
    if isinstance(rgb_val, tuple):
        rgb_hex = mpl.colors.rgb2hex(rgb_val)
    # If multiple values (matrix with rows of rgb values)
    else:
        rgb_hex = []
        for i in range(rgb_val.shape[0]):
            rgb_hex.append(mpl.colors.rgb2hex(rgb_val[i,]))
    # ----------------------------------------------------------------.
    # Set back NaN
    rgb_hex = np.array(rgb_hex)
    if np.sum(idx_nan) > 0:
        rgb_hex[idx_nan] = "NaN"
    # ----------------------------------------------------------------.
    if nan_color is not None:
        rgb_hex[idx_nan] = nan_color
    #################################################
    # If x is dictionary --> Recreate dictionary ####
    #################################################
    if flag_dict:
        rgb_hex = dict(zip(keys, rgb_hex))
    # ---------------------------------------------------------------------.
    return rgb_hex


#### Test get_colors_from_cmap
# get_colors_from_cmap(2, cmap_name="Spectral", vmin=0, vmax = 8)
# get_colors_from_cmap(np.array([2,np.nan,5]), cmap_name="Green", vmin=0, vmax = 8)
# get_colors_from_cmap(np.array([0,2,3,5]), cmap_name="Wistia", vmin=0, vmax = 5)
