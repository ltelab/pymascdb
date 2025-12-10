#!/usr/bin/env python3
"""
Created on Tue Sep 14 11:46:52 2021.

@author: ghiggi
"""
import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar
from skimage import exposure
from skimage.filters import rank
from skimage.morphology import rectangle

# Performs Gamma Correction on the input image.
# skimage.exposure.adjust_gamma(image[, ...])

# Performs Logarithmic correction on the input image.
# skimage.exposure.adjust_log(image[, gain, inv])

# Performs Sigmoid Correction on the input image
# skimage.exposure.adjust_sigmoid(image[, ...])

# http://www.janeriksolem.net/histogram-equalization-with-python-and.html


#################################################
### Workhorse to compute 2D image descriptors ###
#################################################
def _compute_2Dimage_descriptors(da, fun, labels, x="x", y="y", fun_kwargs=None, dask="parallelized"):
    # Checks arguments
    if fun_kwargs is None:
        fun_kwargs = {}
    if not isinstance(da, xr.DataArray):
        raise TypeError("Expecting a xr.DataArray.")
    if not isinstance(x, str):
        raise TypeError("'x' must be a string indicating the width dimension name of the DataArray.")
    if not isinstance(y, str):
        raise TypeError("'y' must be a string indicating the height dimension name of the DataArray.")
    if not isinstance(labels, (str, list)):
        raise TypeError("Descriptor 'labels' must be a str or list of strings.")
    if isinstance(labels, str):
        labels = [labels]
    labels = np.array(labels)
    if not isinstance(labels[0].item(), str):
        raise ValueError("Descriptor 'labels' must be a list of strings.")
    # -----------------------------------------------------------------------.
    # Retrieve DataArray dimension original order
    dims = da.dims
    # -----------------------------------------------------------------------.
    # Check x and y are dimension of the DataArray
    if x not in dims:
        raise ValueError(f"x={x!r} is not a dimension of the DataArray")
    if y not in dims:
        raise ValueError(f"y={y!r} is not a dimension of the DataArray")
    # -----------------------------------------------------------------------.
    ### Retrieve dimensions to eventually stack along a new third dimension
    unstacked_dims = list(set(dims).difference([x, y]))
    # - If only x and y, do nothing
    if len(unstacked_dims) == 0:
        # raise ValueError("Expecting a DataArray with a third dimension in "
        #                  " addition to {!r} and  {!r}".format(x,y))
        stack_dict = {}
        da_stacked = da
    # - If there is already a third dimension, transpose to the last
    elif len(unstacked_dims) == 1:
        img_id = unstacked_dims[0]
        stack_dict = {}
        da_stacked = da.stack(stack_dict).transpose(..., img_id)
    # - If there is more than 3 dimensions, stack it all into a new third dimension
    elif len(unstacked_dims) > 1:
        img_id = "img_id"
        stack_dict = {img_id: unstacked_dims}
        # Stack all additional dimensions into a 3D array with all img_id in the last dimension
        da_stacked = da.stack(stack_dict).transpose(..., img_id)
    else:
        raise NotImplementedError
    # -----------------------------------------------------------------------.
    ### Check the function
    # TODO: checks that len(labels) = len(arr)

    # -----------------------------------------------------------------------.
    ### Compute descriptors for each 2D image
    vectorize = True  # because the function work only on 2D image
    da_stacked = xr.apply_ufunc(
        fun,
        da_stacked,
        input_core_dims=[[x, y]],
        output_core_dims=[["descriptor"]],  # returned data has one dimension
        kwargs=fun_kwargs,
        dask=dask,
        vectorize=vectorize,
        dask_gufunc_kwargs={"output_sizes": {"descriptor": len(labels)}},
        output_dtypes=["float64"],
    )  # TODO: automate

    # Compute the descriptors
    with ProgressBar():
        da_stacked = da_stacked.compute()

    # Add descriptor coordinates
    da_stacked = da_stacked.assign_coords({"descriptor": labels})

    # -----------------------------------------------------------------------.
    # Retrieve dataarray of descriptors
    da_descriptors = da_stacked.unstack(stack_dict)

    # -----------------------------------------------------------------------.
    return da_descriptors


##----------------------------------------------------------------------------.
#####################################
### Workhorse to modify 2D images ###
#####################################
def apply_2Dimage_fun(da, fun, x="x", y="y", fun_kwargs=None):
    """Apply a function to each 2D image in a DataArray.

    This function applies a user-defined function to 2D images stored in a xarray DataArray.
    It handles DataArrays with multiple dimensions by stacking and unstacking as needed,
    ensuring the function is applied to each 2D image independently.

    Parameters
    ----------
    da : xarray.DataArray
        Input DataArray containing 2D images.
    fun : callable
        Function to apply to each 2D image. Should accept a 2D numpy array and return a 2D numpy array.
    x : str, optional
        Name of the width dimension. Default is "x".
    y : str, optional
        Name of the height dimension. Default is "y".
    fun_kwargs : dict, optional
        Additional keyword arguments to pass to the function. Default is None.

    Returns
    -------
    xarray.DataArray
        DataArray with the function applied to each 2D image, maintaining original dimensions.

    Raises
    ------
    TypeError
        If da is not a xarray.DataArray or if x/y are not strings.
    ValueError
        If x or y are not dimensions of the DataArray.
    """
    # Checks arguments
    if fun_kwargs is None:
        fun_kwargs = {}
    if not isinstance(da, xr.DataArray):
        raise TypeError("Expecting a xr.DataArray.")
    if not isinstance(x, str):
        raise TypeError("'x' must be a string indicating the width dimension name of the DataArray.")
    if not isinstance(y, str):
        raise TypeError("'y' must be a string indicating the height dimension name of the DataArray.")
    # -----------------------------------------------------------------------.
    # Retrieve DataArray dimension original order
    dims = da.dims
    # -----------------------------------------------------------------------.
    # Check x and y are dimension of the DataArray
    if x not in dims:
        raise ValueError(f"x={x!r} is not a dimension of the DataArray")
    if y not in dims:
        raise ValueError(f"y={y!r} is not a dimension of the DataArray")
    # -----------------------------------------------------------------------.
    ### Retrieve dimensions to eventually stack along a new third dimension
    unstacked_dims = list(set(dims).difference([x, y]))
    # - If only x and y, do nothing
    if len(unstacked_dims) == 0:
        # raise ValueError("Expecting a DataArray with a third dimension in "
        #                  " addition to {!r} and  {!r}".format(x,y))
        stack_dict = {}
        da_stacked = da
    # - If there is already a third dimension, transpose to the last
    elif len(unstacked_dims) == 1:
        img_id = unstacked_dims[0]
        stack_dict = {}
        da_stacked = da.stack(stack_dict).transpose(..., img_id)
    # - If there is more than 3 dimensions, stack it all into a new third dimension
    elif len(unstacked_dims) > 1:
        img_id = "img_id"
        stack_dict = {img_id: unstacked_dims}
        # Stack all additional dimensions into a 3D array with all img_id in the last dimension
        da_stacked = da.stack(stack_dict).transpose(..., img_id)
    else:
        raise NotImplementedError
    # -----------------------------------------------------------------------.
    ### Apply the function to each 2D image
    dask = "parallelized"  # 'allowed'
    vectorize = True  # because the function work only on 2D image
    da_stacked = xr.apply_ufunc(
        fun,
        da_stacked,
        input_core_dims=[[x, y]],
        output_core_dims=[[x, y]],
        kwargs=fun_kwargs,
        dask=dask,
        vectorize=vectorize,
        output_dtypes=da_stacked.values.dtype,
    )
    # -----------------------------------------------------------------------.
    # Unstack back to original dimensions
    da = da_stacked.unstack(stack_dict).transpose(*dims)

    # -----------------------------------------------------------------------.
    return da


##----------------------------------------------------------------------------.
##################
### Zoom utils ###
##################
def _internal_bbox(img):
    rows = np.any(img, axis=1)
    cols = np.any(img, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    return rmin, rmax, cmin, cmax


def _get_zoomed_image(img):
    rmin, rmax, cmin, cmax = _internal_bbox(img)
    zoom_img = img[rmin : rmax + 1, cmin : cmax + 1]
    return zoom_img


def _center_image(img, nrow, ncol):
    r, c = img.shape
    col_incr = int((ncol - c) / 2)
    row_incr = int((nrow - r) / 2)
    arr = np.zeros((nrow, ncol))
    arr[slice(row_incr, row_incr + r), slice(col_incr, col_incr + c)] = img
    return arr


def xri_zoom(da, x="x", y="y", squared=False):
    """Zoom into 2D images by cropping to non-zero regions and centering.

    This function removes zero-valued borders from images, crops to the smallest bounding box
    containing all non-zero pixels, and centers the result. Optionally creates square images.

    Parameters
    ----------
    da : xarray.DataArray
        Input DataArray containing 2D images.
    x : str, optional
        Name of the width dimension. Default is "x".
    y : str, optional
        Name of the height dimension. Default is "y".
    squared : bool, optional
        If True, output images will be square (same height and width).
        If False, output images maintain their aspect ratio. Default is False.

    Returns
    -------
    xarray.DataArray
        DataArray with zoomed and centered images.

    Raises
    ------
    TypeError
        If da is not a xarray.DataArray or if x/y are not strings.
    ValueError
        If x or y are not dimensions of the DataArray.
    """
    # Checks arguments
    if not isinstance(da, xr.DataArray):
        raise TypeError("Expecting a xr.DataArray.")
    if not isinstance(x, str):
        raise TypeError("'x' must be a string indicating the width dimension name of the DataArray.")
    if not isinstance(y, str):
        raise TypeError("'y' must be a string indicating the height dimension name of the DataArray.")
    # -----------------------------------------------------------------------.
    # Retrieve DataArray dimension original order
    dims = da.dims
    # -----------------------------------------------------------------------.
    # Check x and y are dimension of the DataArray
    if x not in dims:
        raise ValueError(f"x={x!r} is not a dimension of the DataArray")
    if y not in dims:
        raise ValueError(f"y={y!r} is not a dimension of the DataArray")
    # -----------------------------------------------------------------------.
    ### Retrieve dimensions to eventually stack along a new third dimension
    unstacked_dims = list(set(dims).difference([x, y]))

    ## Enforce (..., y, x) dimension order
    da = da.transpose(..., y, x)

    # - If only x and y, do nothing
    if len(unstacked_dims) == 0:
        # raise ValueError("Expecting a DataArray with a third dimension in "
        #                  " addition to {!r} and  {!r}".format(x,y))
        stack_dict = {}
        da_stacked = da
        n_imgs = 0

    # - If there is already a third dimension, transpose to the last
    elif len(unstacked_dims) == 1:
        img_id = unstacked_dims[0]
        stack_dict = {}
        da_stacked = da.stack(stack_dict).transpose(..., img_id)
        # Retrieve number of images
        n_imgs = da_stacked.shape[2]

    # - If there is more than 3 dimensions, stack it all into a new third dimension
    elif len(unstacked_dims) > 1:
        img_id = "img_id"
        stack_dict = {img_id: unstacked_dims}
        # Stack all additional dimensions into a 3D array with all img_id in the last dimension
        da_stacked = da.stack(stack_dict).transpose(..., img_id)
        # Retrieve number of images
        n_imgs = da_stacked.shape[2]
    else:
        raise NotImplementedError
    # -----------------------------------------------------------------------.
    # Extract the list of images
    l_imgs = [da_stacked.values] if n_imgs == 0 else [da_stacked.isel({img_id: i}).values for i in range(n_imgs)]
    # -----------------------------------------------------------------------.
    # Zoom a list of image
    l_zoomed = [_get_zoomed_image(img) for img in l_imgs]

    # Retrieve shape of all zoomed images
    l_shapes = [img.shape for img in l_zoomed]

    # Get number of row an columns of the largest zoomed image
    r_max, c_max = (max(n) for n in zip(*l_shapes))

    # Define size of the zoomed image
    if squared:
        r_max = max([r_max, c_max])
        c_max = max([r_max, c_max])

    # Center all images
    l_zoomed = [_center_image(img, nrow=r_max, ncol=c_max) for img in l_zoomed]

    # Assign it to the DataArray (with new x and y dimensions)
    da_stacked = da_stacked.isel(y=slice(0, r_max), x=slice(0, c_max))
    if n_imgs == 0:
        da_stacked.values = l_zoomed[0]
    else:
        da_stacked.values = np.stack(l_zoomed, axis=-1)

    # -----------------------------------------------------------------------.
    # Unstack back to original dimensions
    da = da_stacked.unstack(stack_dict).transpose(*dims)

    # -----------------------------------------------------------------------.
    return da


##----------------------------------------------------------------------------.
############################################
### Functions to enhance image contrast ####
############################################
def _contrast_stretching(img, pmin=1, pmax=99):
    # Get source dtype
    src_dtype = img.dtype
    # Retrieve img mask
    img_mask = img == 0
    # Compute percentiles
    p2, p98 = np.percentile(img, (pmin, pmax))
    # Perform contrast stretching using percentile values as the intensity range
    img = exposure.rescale_intensity(img, in_range=(p2, p98))
    # Perform contrast stretching using min/max of the dtype as the intensity range
    # img = exposure.rescale_intensity(img, dtype=src_dtype)
    # Change dtype
    img = img.astype(src_dtype)
    # Mask regions that were 0 before stretching
    img[img_mask] = 0
    return img


def _hist_equalization(img, adaptive=False, nbins=256, kernel_size=None, clip_limit=0.03):
    """Apply global or adaptive histogram equalization to an image.

    Parameters
    ----------
    img : arr
       Image array
    adaptive : bool, optional
        If False, it employs the classical histogram equalization.
        If True, it employs Contrast Limited Adaptive Histogram Equalization (CLAHE).
        CLAHE uses histograms computed over different tile regions of the image.
        Local details can therefore be enhanced even in regions that are darker or lighter than most of the image.
        The default is False.
    nbins: int, optional
        Number of bins for image histogram. Note: this argument is ignored for integer images,
        for which each integer is its own bin.
    kernel_size: int or array_like, optional
        Argument used by CLAHE.
        Defines the shape of contextual regions used in the algorithm.
        By default, kernel_size is 1/8 of image height by 1/8 of its width.
    clip_limit: float, optional
        Argument used by CLAHE.
        Clipping limit, normalized between 0 and 1 (higher values give more contrast).
        By default clip_limit=0.01.

    Returns
    -------
    img : arr
        Image array after histogram equalization.

    """
    # Get source dtype
    src_dtype = img.dtype
    # Retrieve img mask
    img_mask = img == 0
    # Retrieve mask for hist equalization
    hist_mask = img > 0
    # Perform histogram equalization
    # - The output values are between 0 and 1 !!!
    if not adaptive:
        img = exposure.equalize_hist(img, nbins=nbins, mask=hist_mask)
    else:
        img = exposure.equalize_adapthist(img, clip_limit=clip_limit, kernel_size=kernel_size, nbins=nbins)
    # Rescale to 0-255
    img = img * 255
    # Change dtype
    img = img.astype(src_dtype)
    # Mask regions that were 0 before stretching
    img[img_mask] = 0
    return img


def _local_hist_equalization(img, footprint=None):
    """Equalize an image using local histograms.

    Parameters
    ----------
    img : arr (uint8, uint16)
       Image array
    footprint: array
        The neighborhood expressed as an ndarray of 1 and 0.
        By default it uses a rectangle of size 1/8 of image height and width
        Custom footprints can be easily generated using skimage.morphology functions
        such as <rectangle, disk, square,star, diamond, octagon,...>

    Returns
    -------
    img : arr
        Image array after equalization.

    """
    # Get source dtype
    src_dtype = img.dtype
    # Define footprint is None
    if footprint is None:
        nrows, ncols = img.shape
        footprint = rectangle(int(nrows / 8), int(ncols / 8))
    # Retrieve img mask
    img_mask = img == 0
    # Retrieve mask for hist equalization
    hist_mask = img > 0
    # Perform local equalization
    img = rank.equalize(np.array(img), footprint, mask=hist_mask)
    # Change dtype
    img = img.astype(src_dtype)
    # Mask regions that were 0 before stretching
    img[img_mask] = 0
    return img


### Wrappers
def xri_contrast_stretching(da, x="x", y="y", pmin=2, pmax=98):
    """Apply contrast stretching to 2D images using percentile-based intensity rescaling.

    Contrast stretching improves image contrast by remapping pixel intensities based on
    specified percentiles, expanding the dynamic range of the image.

    Parameters
    ----------
    da : xarray.DataArray
        Input DataArray containing 2D images.
    x : str, optional
        Name of the width dimension. Default is "x".
    y : str, optional
        Name of the height dimension. Default is "y".
    pmin : float, optional
        Lower percentile for intensity remapping. Default is 2.
    pmax : float, optional
        Upper percentile for intensity remapping. Default is 98.

    Returns
    -------
    xarray.DataArray
        DataArray with contrast-stretched images.

    Notes
    -----
    Zero-valued pixels are preserved and not affected by the stretching operation.
    """
    fun_kwargs = {"pmin": pmin, "pmax": pmax}
    da = apply_2Dimage_fun(da=da, fun=_contrast_stretching, x=x, y=y, fun_kwargs=fun_kwargs)
    return da


def xri_hist_equalization(da, x="x", y="y", nbins=256, adaptive=False, kernel_size=None, clip_limit=0.01):
    """Apply global or adaptive histogram equalization to 2D images.

    Histogram equalization enhances image contrast by redistributing pixel intensities
    to approximate a uniform distribution. Adaptive equalization (CLAHE) computes
    histograms over local tile regions for better enhancement of local details.

    Parameters
    ----------
    da : xarray.DataArray
        Input DataArray containing 2D images.
    x : str, optional
        Name of the width dimension. Default is "x".
    y : str, optional
        Name of the height dimension. Default is "y".
    nbins : int, optional
        Number of bins for image histogram. Ignored for integer images where each
        integer is its own bin. Default is 256.
    adaptive : bool, optional
        If False, uses classical histogram equalization.
        If True, uses Contrast Limited Adaptive Histogram Equalization (CLAHE).
        Default is False.
    kernel_size : int or array_like, optional
        Shape of contextual regions used in CLAHE algorithm.
        By default, uses 1/8 of image height by 1/8 of image width.
    clip_limit : float, optional
        Clipping limit for CLAHE, normalized between 0 and 1.
        Higher values give more contrast. Default is 0.01.

    Returns
    -------
    xarray.DataArray
        DataArray with histogram-equalized images.

    Notes
    -----
    Zero-valued pixels are preserved and not affected by the equalization.
    """
    fun_kwargs = {
        "nbins": nbins,
        "adaptive": adaptive,
        "kernel_size": kernel_size,
        "clip_limit": clip_limit,
    }
    da = apply_2Dimage_fun(da=da, fun=_hist_equalization, x=x, y=y, fun_kwargs=fun_kwargs)
    return da


def xri_local_hist_equalization(da, x="x", y="y", footprint=None):
    """Equalize images using local histograms with a specified neighborhood footprint.

    This function performs histogram equalization using local histograms computed over
    a neighborhood defined by the footprint parameter. This allows for better enhancement
    of local details compared to global histogram equalization.

    Parameters
    ----------
    da : xarray.DataArray
        Input DataArray containing 2D images with dtype uint8 or uint16.
    x : str, optional
        Name of the width dimension. Default is "x".
    y : str, optional
        Name of the height dimension. Default is "y".
    footprint : numpy.ndarray, optional
        The neighborhood expressed as an ndarray of 1's and 0's.
        By default, uses a rectangle of size 1/8 of image height and width.
        Custom footprints can be generated using skimage.morphology functions
        (e.g., rectangle, disk, square, star, diamond, octagon).

    Returns
    -------
    xarray.DataArray
        DataArray with locally equalized images.

    Notes
    -----
    Zero-valued pixels are preserved and not affected by the equalization.
    The input images should have dtype uint8 or uint16.
    """
    fun_kwargs = {"footprint": footprint}
    da = apply_2Dimage_fun(da=da, fun=_local_hist_equalization, x=x, y=y, fun_kwargs=fun_kwargs)
    return da
