"""
Focus measure computation module for MASC analysis.

This module provides focus-measure operators used as textural descriptors
for snowflake ROI analysis.

Functions for:
- Focus and sharpness metrics (LAPM, HISE, WAVS, and other operators)

Translated from fmeasure.m (Said Pertuz 2010) and adapted for Python
Last update : June 2026
Author : Baptiste Carmier
"""

import numpy as np
import pywt

from typing import Optional, Tuple
from scipy import ndimage

def fmeasure(image: np.ndarray, measure: str, roi: Optional[Tuple[int, int, int, int]] = None) -> float:
    """
    Measure the relative degree of focus of an image.

    Applies the requested focus operator to the full image or to a cropped ROI.
    MASC typically uses LAPM, HISE, and WAVS; many other operators from the
    original MATLAB implementation are also supported.

    Args:
        image: Grayscale image (uint8 or float array).
        measure: Focus operator name (e.g. 'LAPM', 'HISE', 'WAVS').
        roi: Optional ROI as (x, y, width, height). None processes the whole image.

    Returns:
        Computed focus value as a float. Returns 0.0 for empty or invalid ROIs.

    Notes:
        - Measure names are case-insensitive.
        - Returns 0.0 instead of raising when the cropped ROI is empty.

    Raises:
        ValueError: If measure is not a supported operator name.
    """
    # Crop to ROI if specified
    if roi is not None:
        x, y, width, height = roi
        image = image[y:y+height, x:x+width]
    
    # Check for empty or invalid image after ROI
    if image.size == 0 or image.shape[0] == 0 or image.shape[1] == 0:
        return 0.0
    
    measure = measure.upper()
    WSize = 15  # Size of local window (used by some operators)
    
    # Using an if-elif structure to handle different measures,
    # while using a switch-case in MATLAB. But Python does not have switch-case natively.
    # A switch-case like was implemented in python >3.10 with match-case, but for compatibility we use if-elif.

    if measure == 'ACMO':
        # Absolute Central Moment (Shirvaikar2004)
        def AcMomentum(values):
            M, N = values.shape
            hist, _ = np.histogram(values.ravel(), bins=256, range=(0, 256))
            hist = hist.astype(float) / (M * N)

            mean_val = np.mean(values)
            hist = np.abs(np.arange(256) - 255*mean_val) * hist
            return float(np.sum(hist))
        
        # Convert to uint8 if needed
        if image.dtype != np.uint8:
            image = (image * 255).astype(np.uint8)

        FM = AcMomentum(image)
        return FM


    elif measure == 'BREN':
        # Brenner's (Santos97)
        M, N = image.shape
        DH = np.zeros_like(image, dtype=float)
        DV = np.zeros_like(image, dtype=float)
        
        if M > 2:
            DH[:M-2, :] = np.diff(image.astype(float), n=2, axis=0)
        if N > 2:
            DV[:, :N-2] = np.diff(image.astype(float), n=2, axis=1)
        
        FM = np.maximum(DH, DV)
        FM = FM ** 2
        FM = float(np.mean(FM))
        return FM
    
    elif measure == 'CONT':
        # Image contrast (Nanda2001)
        def contrast_func(values):
            # values is flattened 3x3 window, center is at index 4
            return np.sum(np.abs(values - values[4]))
        
        FM = ndimage.generic_filter(image.astype(float), contrast_func, size=3, mode='nearest')
        FM = float(np.mean(FM))
        return FM
    
    elif measure == 'CURV':
        # Image Curvature (Helmli2001)
        if image.dtype != np.uint8:
            image = (image * 255).astype(np.uint8)
        
        M1 = np.array([[-1, 0, 1],
                       [-1, 0, 1],
                       [-1, 0, 1]], dtype=float)
        
        M2 = np.array([[1, 0, 1],
                       [1, 0, 1],
                       [1, 0, 1]], dtype=float)
        
        P0 = ndimage.convolve(image.astype(float), M1, mode='nearest') / 6
        P1 = ndimage.convolve(image.astype(float), M1.T, mode='nearest') / 6
        P2 = 3 * ndimage.convolve(image.astype(float), M2, mode='nearest') / 10 - ndimage.convolve(image.astype(float), M2.T, mode='nearest') / 5
        P3 = -ndimage.convolve(image.astype(float), M2, mode='nearest') / 5 + 3 * ndimage.convolve(image.astype(float), M2.T, mode='nearest') / 10 ### P3 peut être faux, avec la trasposition??
        
        FM = np.abs(P0) + np.abs(P1) + np.abs(P2) + np.abs(P3)
        FM = float(np.mean(FM))
        return FM
    
    elif measure == 'DCTE':
        # DCT energy ratio (Shen2006)
        from scipy.fftpack import dct
        
        def dct_ratio(values):
            block = values.reshape(8, 8)
            MT = dct(dct(block.T, norm='ortho').T, norm='ortho') ** 2
            numerator = np.sum(MT) - MT[0, 0]

            if MT[0, 0] != 0:
                return numerator / MT[0, 0]
            else:
                return 0
        
        FM = ndimage.generic_filter(image.astype(float), dct_ratio, size=8, mode='nearest')
        FM = float(np.mean(FM))
        return FM
    
    elif measure == 'DCTR':
        # DCT reduced energy ratio (Lee2009)
        from scipy.fftpack import dct
        
        def re_ratio(values):
            block = values.reshape(8, 8)
            M = dct(dct(block.T, norm='ortho').T, norm='ortho')
            numerator = M[0, 1]**2 + M[0, 2]**2 + M[1, 0]**2 + M[1, 1]**2 + M[2, 0]**2

            if M[0, 0] != 0:
                return numerator / M[0, 0]**2
            else:
                return 0    
        
        FM = ndimage.generic_filter(image.astype(float), re_ratio, size=8, mode='nearest')
        FM = float(np.mean(FM))
        return FM
    
    elif measure == 'GDER':
        # Gaussian derivative (Geusebroek2000)
        N = WSize // 2
        sig = N / 2.5
        x, y = np.meshgrid(np.arange(-N, N+1), np.arange(-N, N+1))
        G = np.exp(-(x**2 + y**2) / (2 * sig**2)) / (2 * np.pi * sig)
        
        Gx = -x * G / (sig**2)
        Gx = Gx / np.sum(Gx)
        
        Gy = -y * G / (sig**2)
        Gy = Gy / np.sum(Gy)
        
        Rx = ndimage.convolve(image.astype(float), Gx, mode='nearest')
        Ry = ndimage.convolve(image.astype(float), Gy, mode='nearest')
        
        FM = Rx**2 + Ry**2
        FM = float(np.mean(FM))
        return FM

    elif measure == 'GLVA':
        # Graylevel variance (Krotkov86)
        FM = float(np.std(image))
        return FM

    elif measure == 'GLLV':
        # Graylevel local variance (Pech2000)
        LVar = ndimage.generic_filter(image.astype(float), np.std, size=WSize, mode='nearest') ** 2
        FM = float(np.var(LVar))
        return FM

    elif measure == 'GLVN':
        # Normalized GLV (Santos97)
        mean_val = np.mean(image)
        if mean_val == 0:
            return 0.0

        FM = float(np.var(image) / mean_val)
        return FM

    elif measure == 'GRAE':
        # Energy of gradient (Subbarao92a)
        Ix = np.zeros_like(image, dtype=float)
        Iy = np.zeros_like(image, dtype=float)
        
        Iy[:-1, :] = np.diff(image.astype(float), axis=0)
        Ix[:, :-1] = np.diff(image.astype(float), axis=1)
        
        FM = Ix**2 + Iy**2
        FM = float(np.mean(FM))
        return FM

    elif measure == 'GRAT':
        # Thresholded gradient (Santos97)
        Th = 0  # Threshold
        
        Ix = np.zeros_like(image, dtype=float)
        Iy = np.zeros_like(image, dtype=float)
        
        Iy[:-1, :] = np.diff(image.astype(float), axis=0)
        Ix[:, :-1] = np.diff(image.astype(float), axis=1)
        
        FM = np.maximum(np.abs(Ix), np.abs(Iy))
        FM[FM < Th] = 0
        
        count = np.sum(FM != 0) # equivalent to np.count_nonzero(FM) 
        if count == 0:
            return 0.0
        
        FM = float(np.sum(FM) / count)
        return FM

    elif measure == 'GRAS':
        # Squared gradient (Eskicioglu95)
        Ix = np.diff(image.astype(float), axis=1)
        FM = Ix ** 2
        FM = float(np.mean(FM))
        return FM

    elif measure == 'HELM':
        # Helmli's mean method (Helmli2001)
        mean_kernel = np.ones((WSize, WSize)) / (WSize * WSize)
        U = ndimage.convolve(image.astype(float), mean_kernel, mode='nearest')
        
        # Avoid division by zero
        R1 = np.divide(U, image.astype(float), out=np.ones_like(U), where=image!=0)
        R1[image == 0] = 1
        
        FM = 1.0 / R1
        FM[U > image] = R1[U > image]
        
        FM = float(np.mean(FM))
        return FM   
            
    # Might be also done with skimage.measure shanon_entropy        
    elif measure == 'HISE':
        # Histogram entropy (Krotkov86)
        if image.size == 0:
            return 0.0
        
        hist, _ = np.histogram(image.ravel(), bins=256, range=(0, 256))
        hist = hist[hist > 0]
        
        if len(hist) == 0:
            return 0.0
        
        prob = hist / hist.sum()
        entropy = -np.sum(prob * np.log2(prob))
        FM = float(entropy)
        return FM

    elif measure == 'HISR':
        # Histogram range (Firestone91)
        FM = float(np.max(image) - np.min(image))
        return FM

    # Might be also done with skimage.filters laplace
    elif measure == 'LAPE':
        # Energy of laplacian (Subbarao92a)
        laplacian = np.array([[0, 1, 0],
                              [1, -4, 1],
                              [0, 1, 0]], dtype=float)
        FM = ndimage.convolve(image.astype(float), laplacian, mode='nearest')
        FM = float(np.mean(FM ** 2))
        return FM

    elif measure == 'LAPM':
        # Modified Laplacian (Nayar89)
        M = np.array([-1, 2, -1], dtype=float)
        Lx = ndimage.convolve(image.astype(float), M.reshape(1, 3), mode='nearest')
        Ly = ndimage.convolve(image.astype(float), M.reshape(3, 1), mode='nearest')
        if np.issubdtype(image.dtype, np.integer):
            # MATLAB imfilter preserves the input class: for uint8 images the
            # convolution output is saturated to [0, 255], so negative
            # responses are zeroed BEFORE abs(), and the uint8 addition
            # abs(Lx)+abs(Ly) saturates at 255 too. Emulate that arithmetic,
            # otherwise the result is ~2x the MATLAB value.
            lo, hi = np.iinfo(image.dtype).min, np.iinfo(image.dtype).max
            Lx = np.clip(Lx, lo, hi)
            Ly = np.clip(Ly, lo, hi)
            FM = np.clip(np.abs(Lx) + np.abs(Ly), lo, hi)
        else:
            FM = np.abs(Lx) + np.abs(Ly)
        FM = float(np.mean(FM))
        return FM
    
    elif measure == 'LAPV':
        # Variance of laplacian (Pech2000)
        laplacian = np.array([[0, 1, 0],
                              [1, -4, 1],
                              [0, 1, 0]], dtype=float)
        ILAP = ndimage.convolve(image.astype(float), laplacian, mode='nearest')
        FM = float(np.var(ILAP))
        return FM
    
    elif measure == 'LAPD': 
        # Diagonal laplacian (Thelen2009)
        M1 = np.array([-1, 2, -1], dtype=float)
        M2 = np.array([[0, 0, -1],
                       [0, 2, 0],
                       [-1, 0, 0]], dtype=float) / np.sqrt(2)
        
        M3 = np.array([[-1, 0, 0],
                       [0, 2, 0],
                       [0, 0, -1]], dtype=float) / np.sqrt(2)
        
        F1 = ndimage.convolve(image.astype(float), M1.reshape(1, 3), mode='nearest')
        F2 = ndimage.convolve(image.astype(float), M2, mode='nearest')
        F3 = ndimage.convolve(image.astype(float), M3, mode='nearest')
        F4 = ndimage.convolve(image.astype(float), M1.reshape(3, 1), mode='nearest')
        
        FM = np.abs(F1) + np.abs(F2) + np.abs(F3) + np.abs(F4)
        FM = float(np.mean(FM))
        return FM
    
    elif measure == 'SFIL':
        # Steerable filters (Minhas2009)
        N = WSize // 2
        sig = N / 2.5
        x, y = np.meshgrid(np.arange(-N, N+1), np.arange(-N, N+1))
        G = np.exp(-(x**2 + y**2) / (2 * sig**2)) / (2 * np.pi * sig)
        
        Gx = -x * G / (sig**2)
        Gx = Gx / np.sum(Gx)
        
        Gy = -y * G / (sig**2)
        Gy = Gy / np.sum(Gy)
    
        
        # Compute responses at different angles
        angles = [0, 45, 90, 135, 180, 225, 270, 315]
        R = np.zeros((image.shape[0], image.shape[1], len(angles)))
        # Weird, missing angle 0 in the original code?
        R[:,:,0] = ndimage.convolve(image.astype(float), Gx, mode='nearest')
        R[:,:,1] = ndimage.convolve(image.astype(float), Gy, mode='nearest')
        R[:,:,2] = np.cos(np.radians(45)) * R[:,:,0] + np.sin(np.radians(45)) * R[:,:,1] ## np.cos() and np.sin() take radians as input
        # Weird, missing angle 90 in the original code?
        R[:,:,3] = np.cos(np.radians(135)) * R[:,:,0] + np.sin(np.radians(135)) * R[:,:,1]
        R[:,:,4] = np.cos(np.radians(180)) * R[:,:,0] + np.sin(np.radians(180)) * R[:,:,1]
        R[:,:,5] = np.cos(np.radians(225)) * R[:,:,0] + np.sin(np.radians(225)) * R[:,:,1]
        R[:,:,6] = np.cos(np.radians(270)) * R[:,:,0] + np.sin(np.radians(270)) * R[:,:,1]
        R[:,:,7] = np.cos(np.radians(315)) * R[:,:,0] + np.sin(np.radians(315)) * R[:,:,1]

        FM = np.max(R, axis=2)
        FM = float(np.mean(FM))
        return FM
    
    elif measure == 'SFRQ':
        # Spatial frequency (Eskicioglu95)
        Ix = np.zeros_like(image, dtype=float)
        Iy = np.zeros_like(image, dtype=float)
        
        Ix[:, :-1] = np.diff(image.astype(float), axis=1)
        Iy[:-1, :] = np.diff(image.astype(float), axis=0)
        
        FM = np.sqrt(Iy**2 + Ix**2)
        FM = float(np.mean(FM))
        return FM

    elif measure == 'TENG':
        # Tenengrad (Krotkov86)
        # Sobel operator
        Sx = np.array([[-1, 0, 1],
                       [-2, 0, 2],
                       [-1, 0, 1]], dtype=float)
        
        Gx = ndimage.convolve(image.astype(float), Sx, mode='nearest')
        Gy = ndimage.convolve(image.astype(float), Sx.T, mode='nearest')
        
        G = Gx**2 + Gy**2
        FM = float(np.mean(G))
        return FM

    elif measure == 'TENV':
        # Tenengrad variance (Pech2000)
        # Sobel operator
        Sx = np.array([[-1, 0, 1],
                       [-2, 0, 2],
                       [-1, 0, 1]], dtype=float)
        
        Gx = ndimage.convolve(image.astype(float), Sx, mode='nearest')
        Gy = ndimage.convolve(image.astype(float), Sx.T, mode='nearest')
        
        G = Gx**2 + Gy**2
        FM = float(np.var(G))
        return FM

    elif measure == 'VOLA':
        # Vollath's correlation (Santos97)
        image_f = image.astype(float)
        
        I1 = np.zeros_like(image_f)
        I2 = np.zeros_like(image_f)
        
        I1[:-1, :] = image_f[1:, :]
        I2[:-2, :] = image_f[2:, :]
        
        FM = image_f * (I1 - I2)
        FM = float(np.mean(FM))
        return FM

    elif measure == 'WAVS':
        # Sum of Wavelet coefficients (Yang2003)
        # Match MATLAB: wavedec2(...,1,'db6') + wrcoef2 (same size as Image).
        # pywt.idwt2 can return +1 row/col vs the input; crop to Image size.
        image_f = image.astype(float)
        h0, w0 = image_f.shape[:2]
        coeffs = pywt.dwt2(image_f, 'db6')
        cA, (cH, cV, cD) = coeffs

        H = pywt.idwt2((None, (cH, None, None)), 'db6')
        V = pywt.idwt2((None, (None, cV, None)), 'db6')
        D = pywt.idwt2((None, (None, None, cD)), 'db6')

        H = H[:h0, :w0]
        V = V[:h0, :w0]
        D = D[:h0, :w0]

        FM = np.abs(H) + np.abs(V) + np.abs(D)
        FM = float(np.mean(FM))
        return FM

    elif measure == 'WAVV':
        # Variance of wavelet coefficients (Yang2003)
        image_f = image.astype(float)
        h0, w0 = image_f.shape[:2]
        coeffs = pywt.dwt2(image_f, 'db6')
        cA, (cH, cV, cD) = coeffs

        H = np.abs(pywt.idwt2((None, (cH, None, None)), 'db6'))[:h0, :w0]
        V = np.abs(pywt.idwt2((None, (None, cV, None)), 'db6'))[:h0, :w0]
        D = np.abs(pywt.idwt2((None, (None, None, cD)), 'db6'))[:h0, :w0]

        FM = float(np.var(H) + np.var(V) + np.var(D))
        return FM
    
    elif measure == 'WAVR':
        # Wavelet ratio - ratio of high-pass to low-pass energies
        coeffs = pywt.wavedec2(image.astype(float), 'db6', level=3) # use of wavedec2 to get multiple levels, as level 3 is mentionned in matlab original code
        
        # coeffs structure: [cA_n, (cH_n, cV_n, cD_n), ..., (cH_1, cV_1, cD_1)]
        # coeffs[0] = approximation at level 3 (array)
        # coeffs[1] = (cH3, cV3, cD3) details at level 3
        # coeffs[2] = (cH2, cV2, cD2) details at level 2
        # coeffs[3] = (cH1, cV1, cD1) details at level 1

        # Reconstruct approximation at level 3
        cA3 = coeffs[0]
        A3 = np.abs(coeffs[0])  # an array, not tuple

        # Reconstruct approximation at level 2
        cH3, cV3, cD3 = coeffs[1]
        A2 = np.abs(pywt.wrcoef2((cA3, (cH3, cV3, cD3)), 'db6'))

        # Reconstruct approximation at level 1
        cH2, cV2, cD2 = coeffs[2] 
        A1 = np.abs(pywt.waverec2((A2, (cH2, cV2, cD2)), 'db6'))

        # Reconstruct details at level 1 
        cH1, cV1, cD1 = coeffs[3]
        H1 = np.abs(pywt.idwt2((None, (cH1, None, None)), 'db6'))
        V1 = np.abs(pywt.idwt2((None, (None, cV1, None)), 'db6'))
        D1 = np.abs(pywt.idwt2((None, (None, None, cD1)), 'db6'))

        print(A1.shape, A2.shape, A3.shape)
        A = A1 + A2 + A3

        # Calculate high-frequency energy from finest detail coefficients
        WH = H1**2 + V1**2 + D1**2
        WH = np.mean(WH)
        WL = np.mean(A)
        
        if WL == 0:
           return 0.0
        
        FM = float(WH / WL)
        return FM

    else:
        raise ValueError(f"Unknown measure: {measure}")