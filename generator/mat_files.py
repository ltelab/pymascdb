import numpy as np
from scipy import io as sio
from datetime import datetime, timedelta 
import numpy as np


import pandas as pd
import matplotlib.pyplot as plt 

def digits_dictionary():
    """
    Dictionary containing the number of decimal places to be used for various variables
    """
    dict= {  'index':         0,
             'pix_size':      7,
             'quality_xhi':   1,

             'n_roi'   :      0,
             'flake_n_roi'  : 0,
             'area'    :      9,
             'perim'   :      5,
             'Dmean'   :      5,
             'Dmax'    :      5,   
             'eq_radius':     5,
             'area_porous':   9,
             'area_porous_r': 3,
             'ell_fit_A':     5,
             'ell_fit_B':     5,
             'ell_fit_area':  9,
             'ell_fit_ori':   1,
             'ell_fit_ecc':   2,
             'compactness':   2,
             'ell_in_A':      5,
             'ell_in_B':      5,
             'ell_in_area':   9,
             'ell_out_A':     5,
             'ell_out_B':     5,
             'ell_out_area':  9,
             'roundness':     2,
             'p_circ_out_r':  2,
             'rectangularity':2,
             'bbox_width':    5,
             'bbox_len':      5,
             'rect_perim_ratio':  2,
             'rect_aspect_ratio': 2,
             'rect_eccentricity': 2,
             'solidity':          2,
             'convexity':         2,
             'hull_n_angles':     0,
             'p_circ_r':          2,
             'frac_dim_boxcounting': 2,
             'frac_dim_theoretical': 2,
             'nb_holes':             0,
             'skel_N_ends':       0,
             'skel_N_junc' :      0,
             'skel_perim_ratio':  2,
             'skel_area_ratio':   3,
             'sym_P1':            2,
             'sym_P2':            2,
             'sym_P3':            2,
             'sym_P4':            2,
             'sym_P5':            2,
             'sym_P6':            2,
             'sym_Pmax_id':       0,
             'sym_P6_max_ratio':  2,
             'sym_mean':          1,
             'sym_std':           1,
             'sym_std_mean_ratio': 2,
             'intensity_mean':     2,
             'intensity_max':      2,
             'contrast':           2,
             'intensity_std':      1,
             'hist_entropy':       1,
             'local_std':          1,
             'local_intens':       2,
             'lap_energy':         1,
             'wavs':               1,
             'complexity':         2,
             'har_energy':         6,
             'har_contrast':       1,
             'har_corr':           2,
             'har_hom':            3,
             'roi_centroid_X':     0,
             'roi_centroid_Y':     0,
             'roi_width':          0,
             'roi_height':         0,
             'Dmax_ori':           1,
             'Dmax_90':            5,
             'D90_r':              2,
             'riming_class_prob':  2,
             'riming_deg_level':   2,
             'melting_class_id':   0,
             'melting_prob':       2,
             'snowflake_class_prob': 2,

             'fallspeed':          3,
             'latitude':           4,
             'longitude':          4,
             'altitude':           1,
             'flake_quality_xhi':  1,
             'flake_Dmax':         5,

             'gan3d_mass':         9,
             'gan3d_volume':       15,
             'gan3d_gyration':     5,

             'bs_normalized_angle':   3,
             'bs_mixing_ind':         2,

             'env_T':              1,
             'env_P':              1,
             'env_DD':             1,
             'env_FF':             1,
             'env_RH':             1
    }

    return dict


def lat_lon_alt(campaign):
    """
    Get lat, lon (WGS84) and alt. for a given campaign 
    """

    if campaign == 'APRES3-2016':
        lat = -66.6628
        lon = 140.0014
        alt = 41.0 
    elif campaign == 'APRES3-2017':
        lat = -66.6628
        lon = 140.0014
        alt = 41.0 
    elif campaign == 'Davos-2015':
        lat = 46.8297
        lon = 9.8093
        alt = 2540. 
    elif campaign == 'Davos-2019':
        lat = 46.8450
        lon = 9.8716
        alt = 1512.0
    elif campaign == 'ICEGENESIS-2021':
        lat =  47.0830
        lon =  6.7922
        alt = 1018.0 
    elif campaign == 'ICEPOP-BKC-2018': #BoKwang1–ri Community center (before 21.02.2018)
        lat = 37.7382
        lon = 128.7586
        alt = 175.
    elif campaign == 'ICEPOP-2018': #Mayhills supersite if not specified otherwise
        lat = 37.6652
        lon = 128.6996
        alt = 789.
    elif campaign == 'Jura-2019': # Les Charbonnieres
        lat =  46.6702 
        lon =  6.3125
        alt =  1045.0 
    elif campaign == 'Jura-2023': # Les Charbonnieres
        lat =  46.6702 
        lon =  6.3125
        alt =  1045.0 
    elif campaign == 'PLATO-2019':
        lat = -68.5752 
        lon = 77.9659 
        alt = 10.0
    elif campaign == 'POPE-2020':
        lat = -71.9499 
        lon =  23.3471
        alt =  1382.0
    elif campaign == 'Valais-2016':
        lat = 46.1222 
        lon = 7.2122 
        alt = 2370.0  
    elif campaign == 'Remoray-2022':
        lat = 46.7665 
        lon = 6.2402
        alt = 920.0
    elif campaign == 'Norway-2016':
        lat = 59.8118
        lon = 7.2143
        alt = 991.0
    else:
        print("Warning: campaign "+campaign+" not recognized")
        lat = np.nan
        lon = np.nan
        alt = np.nan

    return lat,lon,alt


def datenum_to_datetime(datenum):
    """
    Convert Matlab datenum into Python datetime.
    :param datenum: Date in datenum format
    :return:        Datetime object corresponding to datenum.
    """
    days = datenum % 1
    return datetime.fromordinal(int(datenum)) + timedelta(days=days) - timedelta(days=366)

def pmax06(mat_in):

    vec=[mat_in['P0'][0][0],
         mat_in['P1'][0][0],
         mat_in['P2'][0][0],
         mat_in['P3'][0][0],
         mat_in['P4'][0][0],
         mat_in['P5'][0][0],
         mat_in['P6'][0][0]]
            
    return [np.argmax(vec),np.max(vec)]

def compute_riming_id(prob_vec,use_all_probs=False):
    if use_all_probs:
        return np.sum([1,2,3,4,5]*prob_vec)
    else:
        return np.argmax(prob_vec)+1

def compute_riming_idx(R):
    return 0.5*(np.sin(0.25*np.pi*(R-3))+1)

def get_riming_name(R):
    if R == 1:
        return 'unrimed'
    elif R == 2:
        return 'rimed'
    elif R == 3:
        return 'densely_rimed'
    elif R == 4:
        return 'graupel-like'
    elif R == 5:
        return 'graupel'

def get_melting_name(M):
    if M == 0:
        return 'dry'
    elif M == 1:
        return 'melting'

def pad_with_zeros(A, r=1024, c=1024):
   """
   Roughly center the image A into a r x c grid
   """
   out = np.zeros((r, c),dtype=np.uint8)
   r_, c_ = np.shape(A)

   ll = np.int((r-r_)*0.5)
   bb = np.int((c-c_)*0.5) 

   if (ll >= 0) and (bb >= 0):
        out[ll:(ll+r_), bb:(bb+c_)] = A
        return out
   elif (ll < 0) and (bb > 0): # crop on the right (arbitrary choice)
        out[:,bb:(bb+c_)] = A[0:r,:]
        return out
   elif (ll > 0):              # crop on the top
        out[ll:(ll+r_),:] = A[:,0:c]
        return out
   else:
        out[ll:(ll+r_),bb:(bb+c_)] = A[0:r,0:c]
        return out

def id2name(id_in):

    if id_in == 1:
        return 'small_particle'
    elif id_in == 2:
        return 'columnar_crystal'
    elif id_in == 3:
        return 'planar_crystal'
    elif id_in == 4:
        return 'aggregate'
    elif id_in == 5:
        return 'graupel'
    elif id_in == 6:
        return 'columnar_planar_combination'

def masc_mat_file_to_dict(fn,pix_size=33.5e-6):

    """
    Converts the content of a .mat file obtained with the method
    of Praz et al 2017
    """

    mat=(sio.loadmat(fn))['roi'][0,0]

    # Create a large dictionary to store the .mat data originally created
    # with the Matlab codes of Prazet al 2017
    #   Not all the data are kept as
    # (1) Many descriptors can be computed by combination of others
    # (2) Not all the descriptors are really informative

    # [pix] indicates the 1D pixel size

    #----------------------------------------------------

    # Variables to convert pixels to m
    p1 = pix_size   # to convert [pix] to m
    p2 = p1**2      # to convert [pix**2] to m**2

    #Variables to compute only once
    riming_id=compute_riming_id(mat['riming_probs'][0]) 
    riming_deg_level = compute_riming_idx(compute_riming_id(mat['riming_probs'][0],use_all_probs=True))

    dict={
        # Time
        'datetime': datenum_to_datetime(mat['tnum'][0][0]), #datetime obj
        'flake_id': (fn.split('/')[-1]).split('_cam')[0],
        'flake_number_tmp':   mat['id'][0][0],                          # Temporary flake number (reset to 1 after reboot)
        'pix_size': pix_size,                                           # Pixel size in m

        # Other MASC features   
        'n_roi':                mat['n_roi'][0][0],                     # Number of particles on images, including the main ROI of the feat. above [-] 
        'cam_id':               mat['cam'][0][0],                       # camera id (0,1,2)   


        ## Features from ROI of Praz et al 2017. Numbered as in Table A1.
        # Labeled as I, II, III, -
        # I   = used in hydro classification
        # II  = used in riming degree estimation
        # III = used in melting estimation
        # -   = discarded

        # C1: particle size and area
        'area':         mat['area'][0][0]*p2,          # ROI area  [m**2]                       1. - 
        'perim':        mat['perim'][0][0]*p1,         # ROI perimeter [m]                      2. I, II
        'Dmean':        mat['Dmean'][0][0]*p1,         # Dmean [m] (ROI 0.5*(width x height))   3. -
        'Dmax':         mat['Dmax'][0][0]*p1,          # Dmax  [m]                              4. I
        'eq_radius':    mat['eq_radius'][0][0]*p1,     # Eq. area radius [m]                    5. -
        'area_porous':  mat['area_porous'][0][0]*p2,   # ROI Area with holes removed [m**2]     6. II, III
        'area_porous_r':mat['area_porous'][0][0]/mat['area'][0][0],              #              7. III

        # C2 elliptical approximation
        'ell_fit_A':  mat['E']['a'][0,0][0,0]*p1,         # Fitted ell maj dim [m]                                8. II
        'ell_fit_B':  mat['E']['b'][0,0][0,0]*p1,         # Fitted ell min dim [m]                                9. I,III
        'ell_fit_area':   mat['E']['a'][0,0][0,0]*mat['E']['b'][0,0][0,0]*np.pi*p2,  #[m**2]                      10. -
        'ell_fit_ori':    mat['E']['theta'][0,0][0,0],    # Fitted ell orientation [°]                            11. -
        'ell_fit_a_r':    mat['E']['b'][0,0][0,0]/mat['E']['a'][0,0][0,0], # Fitted ell asp. ratio                12. III
        'ell_fit_ecc':    
                    np.sqrt(1-mat['E']['b'][0,0][0,0]/mat['E']['a'][0,0][0,0]), # Fitted ell eccentric.          13. I
        'compactness':    mat['compactness'][0][0],    # (?) Proj, area /fitted ellipses area ratio [-]          14. I,III
        'ell_in_A':  mat['E_in']['a'][0,0][0,0]*p1,       # Inner ell maj dim [m] --same center as fitted one    15. - 
        'ell_in_B':  mat['E_in']['b'][0,0][0,0]*p1,       # Inner ell min dim [m] --same center as fitted one    16. -
        'ell_in_area':   mat['E_in']['a'][0,0][0,0]*mat['E_in']['b'][0,0][0,0]*np.pi*p2,  #[m**2]                17. I
        'ell_out_A':  
        mat['E_out']['a'][0,0][0,0]*p1,     # Outer ell maj dim [m] --same center as fitted one, orientation as inner one  18. -
        'ell_out_B':  
        mat['E_out']['b'][0,0][0,0]*p1,     # Outer ell min dim [m] --same center as fitted one, orientation as inner      19. I, II, III
        'ell_out_area':   mat['E_out']['a'][0,0][0,0]*mat['E_out']['b'][0,0][0,0]*np.pi*p2,  #[m**2]                       20. -
        
        # C3: particle shape
        'roundness':            mat['roundness'][0][0],                   # area / circum. circle ratio [-]                     30. III
        'p_circ_out_r':         mat['perim'][0][0]/(2.*np.pi*mat['C_out']['r'][0,0][0,0]), # perim /circum p ratio [-]          31. II, III
        'rectangularity':       mat['Rect']['A_ratio'][0,0][0,0],         # Area to bounding box area ratio [-]                 32. -
        'bbox_width':           mat['Rect']['width'][0,0][0,0]*p1,        # Bound. box width [m]                                33. -
        'bbox_len':             mat['Rect']['length'][0,0][0,0]*p1,       # Bound. box length [m]                               34. -
        'rect_perim_ratio':     mat['Rect']['p_ratio'][0,0][0,0],       # Box perimeter to particule perim ratio [-]            35. I, II, III
        'rect_aspect_ratio':    mat['Rect']['aspect_ratio'][0,0][0,0],  # Bounding box aspect ratio [-]                         36. I
        'rect_eccentricity':    mat['Rect']['eccentricity'][0,0][0,0],  # Bounding bo eccentricity [-]                          37. -
        'solidity':             mat['hull']['solidity'][0,0][0,0],      # area / CH area ratio    [-]                           38. I
        'convexity':            mat['hull']['convexity'][0,0][0,0],     # perim / CH perim ratio  [-]                           39. I, III
        'hull_n_angles':        mat['hull']['xh'][0,0].shape[0],    # number of CH vertices   [-]                           40. II, III
        'p_circ_r':             mat['perim'][0][0]/(2.*np.pi*mat['eq_radius'][0][0]), # perim /eq. circ p ratio [-]         41. III
        'frac_dim_boxcounting': mat['F'][0][0],                     # Boxcouting fractal dim  [-]                           42. -
        'frac_dim_theoretical': mat['F_jac'][0][0],                 # Theoretical fractal dim [-] Grazioli et al 2014       43. II
        'nb_holes':             mat['nb_holes'][0][0],              # Number of holes [-]                                   ---

        # C4: morphological skeleton
        'skel_N_ends':      mat['skel']['N_ends'][0,0][0,0],      # Number of skeleton ends [-]                                 44. III
        'skel_N_junc':      mat['skel']['N_junctions'][0,0][0,0], # Number of skeleton junctions [-]                            45. -
        'skel_perim_ratio': mat['skel']['p_ratio'][0,0][0,0],     # Skeleton length to perimeter ratio [-]                      46. -
        'skel_area_ratio':  mat['skel']['A_ratio'][0,0][0,0],     # Skeleton length to area ratio [pix**-1]                     47. -

        # C5: rotational symmetry
        #'sym_P0':       mat['Sym']['P0'][0,0][0,0],                   # Standardized distance to centroid Fourier power spectrum comp. P0 48. -
        'sym_P1':       mat['Sym']['P1'][0,0][0,0],                   # Standardized distance to centroid Fourier power spectrum comp. P1 49. I
        'sym_P2':       mat['Sym']['P2'][0,0][0,0],                   # Standardized distance to centroid Fourier power spectrum comp. P2 50. II
        'sym_P3':       mat['Sym']['P3'][0,0][0,0],                   # Standardized distance to centroid Fourier power spectrum comp. P3 51. I, II, III
        'sym_P4':       mat['Sym']['P4'][0,0][0,0],                   # Standardized distance to centroid Fourier power spectrum comp. P4 52. II, III
        'sym_P5':       mat['Sym']['P5'][0,0][0,0],                   # Standardized distance to centroid Fourier power spectrum comp. P5 53. -
        'sym_P6':       mat['Sym']['P6'][0,0][0,0],                   # Standardized distance to centroid Fourier power spectrum comp. P6 54. I, II
        'sym_Pmax_id':     pmax06(mat['Sym'])[0],                     # Max ID among P0 to P6                                             55. I, III
        'sym_P6_max_ratio': mat['Sym']['P6'][0,0][0,0]/pmax06(mat['Sym'])[1],#                                                            56. II
        'sym_mean':     mat['Sym']['mean'][0,0][0,0],                 # Mean distance to centroid [pix] ?                                 57. 
        'sym_std':      mat['Sym']['std'][0,0][0,0],                  # STD distance to centroid  [pix] ?                                 58. I, II
        'sym_std_mean_ratio': mat['Sym']['std'][0,0][0,0]/mat['Sym']['mean'][0,0][0,0],                              #                    59. I, II

        # C6: texture operators
        'intensity_mean':       mat['mean_intens'][0][0],            # Mean intensity (not sure about units, seems normalized)        60. I, II, III
        'intensity_max':        mat['max_intens'][0][0],             # Max intensity /brightness                                      61. II, III
        'contrast':             mat['contrast'][0][0],               # contrast                                                       62. I
        'intensity_std':        mat['std'][0][0],                    # std intensity                                                  63. III
        'hist_entropy':         mat['hist_entropy'][0][0],           # Brightness histogram entropy                                   64. -
        'local_std':            mat['local_std'][0][0],              # Average grey-level Local standard deviation 3x3                65. I, II, III
        'local_intens':         mat['range_intens'][0][0],           # Local average range intensity (probably 3x3 window)            66. I, III
        'lap_energy':           mat['lap'][0][0],                    # Energy of the Laplacian                                        67. -
        'wavs':                 mat['wavs'][0][0],                   # Sum of wavelet coeff                                           68. I, II
        'complexity':           mat['complex'][0][0],                # Complexity from Garret and Yuter 2014 [-]                      69. II

        # C7: Haralick features / co-occurrence matrix
        'har_energy':           mat['H']['Energy'][0,0][0,0],             # Haralick Energy      ? 70. -
        'har_contrast':         mat['H']['Contrast'][0,0][0,0],           # Haralick contrast      71. II, III
        'har_corr':             mat['H']['Correlation'][0,0][0,0],        # Haralick correlation   72. I, II, III
        'har_hom':              mat['H']['Homogeneity'][0,0][0,0],        # Haralick homogeneity   73., I, II

        # ------------------------- END of Praz 2017 input features -------------------------------------
        # Other output or direct features

        'roi_centroid_X':           mat['centroid'][0][0],              # Centroid X-pos of ROI -with respect to the raw picture
        'roi_centroid_Y':           mat['centroid'][0][1],              # Centroid Y-pos of ROI     
        'roi_width':                mat['width'][0][0],                 # ROI x size [pix]                        
        'roi_height':               mat['height'][0][0],                # ROI y size [pix]                     

        # Quality features
        'quality_xhi':              mat['xhi'][0][0],                   # Quality index. As in Praz et al 2017

        # Other orientations
        'Dmax_ori':             mat['Dmax_theta'][0][0],            # [°] orientation of Dmax

        # Other characteristicts from cross-max dimension (like Baker Lawson 2006)
        'Dmax_90':              mat['D90']['Dmax_90'][0,0][0,0]*p1,              # Maximum dimension in the hort. with respect to Dmax [m]
        #'Dmax_0_90':            mat['D90']['Dmax_0'][0,0][0,0],                 # Test variable. Dmax recalculated from rotation of Dmax_90
        'D90_r':                mat['D90']['AR'][0,0][0,0],                      # Axis ratio of D90 and D0 

        # Riming probabilities for each class 1-5
        'riming_class_id':       riming_id,                       # 1 to 5 
        'riming_class_prob':  round(mat['riming_probs'][0][riming_id-1],2), 
        'riming_deg_level':      round(riming_deg_level,2),       # 0 to 1
        'riming_class_name':     get_riming_name(riming_id),      # Unrimed, rimed, densely_rimed, graupel-like, graupel

        'melting_class_id':      mat['melting_ID'][0][0],                                            # 0 or 1    
        'melting_prob':          mat['melting_probs'][0][0],                                        # 0 to 1
        'melting_class_name':    get_melting_name(mat['melting_ID'][0][0]),
        
        # Hydrometeor classification 
        # 1 = small particle (SP)
        # 2 = columnar crystal (CC)
        # 3 = planar crystal (PC) 
        # 4 = aggregate (AG)
        # 5 = graupel (GR)
        # 6 = combination of columnar and planar crystals (CPC) 
        
        'snowflake_class_name':           id2name(mat['label_ID'][0][0]),                           # short name of hydro class
        'snowflake_class_id':             mat['label_ID'][0][0],                                    # Label ID 1 to 6
        'snowflake_class_prob':        round(mat['label_probs'][0][mat['label_ID'][0][0]-1],2)}  # Prob of label X

    # If melting then riming is undefined
    if mat['melting_ID'][0][0] == 1:
        dict['riming_class_id']       = 0
        dict['riming_class_name']     = 'undefined'
        dict['riming_class_prob']  =  np.nan
        dict['riming_deg_level']      =  np.nan
        
    return dict

def masc_mat_triplet_to_dict(fnames,pix_size=33.5e-6,campaign=''):
    """
    Get or compute from a triplet of .mat files, generated 
    with the method of Praz et al 2017, the set of descriptors,
    properties, or retrievals that are valid for the triplet as 
    a whole. 
    
    For example: fall speed, hydro class, riming degree, mass/volume
    estimates from 3D-GAN, etc. 

    Note: it assumes that a triplet is available. 

    Input: 

    fnames: array with the fname of the triplet. Must have length 3

    Output: 

    dictionary

    """

    if len(fnames) != 3:
        print("Exactly 3 file names must be provided")
        return None

    # Read the 3 files
    mat0=(sio.loadmat(fnames[0]))['roi'][0,0]
    mat1=(sio.loadmat(fnames[1]))['roi'][0,0]
    mat2=(sio.loadmat(fnames[2]))['roi'][0,0]

    # Check if fallspeed was recorded
    try:
        fs=mat0['fallspeed'][0][0]
    except:
        fs=np.nan

    # Get quality of the images
    Xhi=np.mean([mat0['xhi'][0][0],mat1['xhi'][0][0],mat2['xhi'][0][0]])

    # Get Dmax
    Dmax=np.max([mat0['Dmax'][0][0],mat1['Dmax'][0][0],mat2['Dmax'][0][0]])*pix_size

    # Get probabilities of hydrometeor classif to assign the proper class
    probs      = mat0['label_probs'][0]+mat1['label_probs'][0]+mat2['label_probs'][0]
    label_id   = np.argmax(probs)+1
    label_prob = np.max(probs/np.sum(probs))

    # Get riming degrees
    probs = mat0['riming_probs'][0]+mat1['riming_probs'][0]+mat2['riming_probs'][0]
    probs = probs/np.sum(probs)
    riming_id               = compute_riming_id(probs)                                          # Index 1 to 5
    riming_id_prob          = probs[riming_id-1]
    riming_deg_level        = compute_riming_idx(compute_riming_id(probs,use_all_probs=True))   # Index 0 to 1 

    # Melting
    melting_prob =(mat0['melting_probs'][0]+mat1['melting_probs'][0]+mat2['melting_probs'][0])/3.
    melting_id   =np.round(melting_prob)

    # Lat, Lon, alt
    lat, lon, alt = lat_lon_alt(campaign) 

    dict={
        # Time and location
        'datetime' : datenum_to_datetime(mat0['tnum'][0][0]), # datetime obj
        'campaign' : campaign,
        'latitude' :lat,  # °N WGS84
        'longitude':lon,  # °E WGS84
        'altitude' :alt,  #  m. a. msl

        # Flake info
        'flake_id': (fnames[0].split('/')[-1]).split('_cam')[0],   # Flake ID unique
        'flake_number_tmp':   mat0['id'][0][0],                    # Temporary flake number (reset to 1 after reboot)
        'flake_quality_xhi':      Xhi,                             # Average quality index

        # GLobal info
        'fallspeed': fs,                                     # [m s**-1] fall speed
        'flake_n_roi':  np.mean([mat0['n_roi'][0][0],mat1['n_roi'][0][0],mat2['n_roi'][0][0]]),  # Avg # of particles per cam [-] 

        # Dmax
        'flake_Dmax':    Dmax,                               # m

        # Riming degree
        'riming_deg_level':     round(riming_deg_level,2),   # 1 to 5 
        'riming_class_id':      riming_id,                   # 0 to 1
        'riming_class_prob': round(riming_id_prob,2),
        'riming_class_name':    get_riming_name(riming_id),

        # Mmelting
        'melting_class_id': melting_id[0],                   # 0 or 1
        'melting_prob':     melting_prob[0],
        'melting_class_name':    get_melting_name(melting_id[0]),

        # --

        # Hydrometeor classification 
        # 1 = small particle (SP)
        # 2  = columnar crystal (CC)
        # 3 = planar crystal (PC) 
        # 4 = aggregate (AG)
        # 5 = graupel (GR)
        # 6 = combination of columnar and planar crystals (CPC) 
        'snowflake_class_name':              id2name(label_id),               # label id to label name           
        'snowflake_class_id':                label_id,                        # 1 to 6 
        'snowflake_class_prob':           round(label_prob,2),

        # Placeholder for 3D-GAN products of Leinonen et al 2021
        'gan3d_mass':         np.nan,        # mass from 3d-gan [kg]
        'gan3d_volume':       np.nan,        # convex hull  volume from 3d-gan in [m**3]
        'gan3d_gyration':     np.nan,        # gyration radius  [m]

        # Placeholder for blowing snow detection and classification Schaer et al 2020
        # parameter given is the normalized angle and the mixing index (defined only for mixed BS and precip), both 
        # varying between 0 and 1
        # Mixing index closer to 1 indicates more BS.
        # Normalized angle < 0.193 means precip, > 0.881 means BS

        # Placeholder for environmental information 
        'env_T':          np.nan,    # Temperature       [°C]
        'env_P':          np.nan,    # Pressure          [hPa]
        'env_DD':         np.nan,    # Wind direction    [°]
        'env_FF':         np.nan,    # Wind speed        [m/s]
        'env_RH':         np.nan,     # Relative Humidity [%]

        'bs_normalized_angle': np.nan,           # Normalized angle [-]
        'bs_mixing_ind':   np.nan,               # Mixing index [-] 
        'bs_precip_class_name': 'undefined'      # Precipitation type (blowing_snow, precip, mixed or undefined )

    }

    # If melting then riming is undefined
    if dict['melting_class_id'] == 1:
        dict['riming_class_id']       = 0
        dict['riming_class_name']     = 'undefined'
        dict['riming_class_prob']     =  np.nan
        dict['riming_deg_level']      =  np.nan

    return dict

def triplet_images_reshape(fnames,pix_size=33.5e-6,newshape=[1024,1024]):
    """
    Read the filenames of a processed MASC triplet and center them into a 
    grid of common size. Add also a few information about the triplet 
    itself (datetime, flake ID, pixel size)
    """

    if len(fnames) != 3:
        print("Exactly 3 file names must be provided")
        return None

    # Read the 3 files
    A0=pad_with_zeros((sio.loadmat(fnames[0]))['roi'][0,0]['data'],r=newshape[0],c=newshape[0])
    A1=pad_with_zeros((sio.loadmat(fnames[1]))['roi'][0,0]['data'],r=newshape[0],c=newshape[0])
    A2=pad_with_zeros((sio.loadmat(fnames[2]))['roi'][0,0]['data'],r=newshape[0],c=newshape[0])

    return np.dstack([A0,A1,A2])

    return dict



