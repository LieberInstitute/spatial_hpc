import os
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')

from pyhere import here
from pathlib import Path
import session_info

import numpy as np
import pandas as pd
import json
import sys
from loopy.sample import Sample
import tifffile
from PIL import Image
import re
import matplotlib.pyplot as plt

import scanpy as sc
from rasterio import Affine
from loopy.utils.utils import remove_dupes, Url
import re

Image.MAX_IMAGE_PIXELS = None

#this_donor = "Br8325"
this_donor = "Br3942"

sample_info_path = here('processed-data', 'VSPG_image_stitching', 'transformations.csv')
combined_dir = Path(here('processed-data', 'VSPG_image_stitching',f'combined_{this_donor}'))
combined_dir.mkdir(parents = True, exist_ok = True)

img_out_fullres = Path(here(combined_dir, f'combined_{this_donor}.tif'))
tissue_out_path = Path(here(combined_dir, f'tissue_positions_{this_donor}.csv'))
img_out_lowres = Path(here(combined_dir, f'combined_tissue_lowres_{this_donor}.png'))
json_out_path = Path(here(combined_dir, 'scalefactors_json.json'))

#   55-micrometer diameter for Visium spot but 65-micrometer spot diameter used
#   in 'spot_diameter_fullres' calculation for spaceranger JSON. The
#   difference between 55 and 65 does indeed exist and is properly documented,
#   but is likely a bug in the sense the choice was probably unintentional
#   https://kb.10xgenomics.com/hc/en-us/articles/360035487812-What-is-the-size-of-the-spots-on-the-Visium-Gene-Expression-Slide-
#   https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/spatial

SPOT_DIAMETER_M = 55e-6
SPOT_DIAMETER_JSON_M = 65e-6

LOWRES_MAX_SIZE = 1200
BACKGROUND_COLOR = 0

#   Read in sample info and subset to samples of interest
sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info = sample_info.loc[
    (sample_info['Brain'] == this_donor) &
    ~sample_info['xml_path'].isna(),# &
    #sample_info['In analysis'],
    :
]

################################################################################
#   Functions
################################################################################

#   Given the size (shape) of a PIL Image 'img', simulate a rotation
#   counter-clockwise by theta (radians) about the top left of the image (0, 0),
#   and return the translation (shape (2,) numpy array) needed to add to the
#   indices of the numpy array given by np.array(img.rotate([theta in degrees]))
#   to correctly apply the rotation about the top left of the image
def adjust_trans(img_shape, theta):
    #   (Negative) angle the top edge of the image makes with the diagonal,
    #   and the length of that diagonal. We'll need these to calculate
    #   where the bottom-right corner of the image moves as a function of
    #   theta
    theta_r = np.arctan(img_shape[1] / img_shape[0])
    r = np.linalg.norm(img_shape)
    #   Consider the image as a rectangle, and rotate that rectangle about the
    #   top-left corner (which is considered the origin). Track where all 4
    #   corners of the rectangle go, since that determines where the originally
    #   top-left corner ends up relative to the expanded-dimension numpy array's
    #   top-left corner. 'adjust' is simply the vector difference between those
    #   corners
    adjust = np.array(
        [
            -1 * max(
                0,
                img_shape[0] * np.sin(theta),
                img_shape[1] * np.sin(theta - np.pi / 2),
                r * np.sin(theta - theta_r)
            ),
            min(
                0,
                img_shape[0] * np.cos(theta),
                img_shape[1] * np.cos(theta - np.pi / 2),
                r * np.cos(theta - theta_r)
            )
        ]
    )
    return adjust

#   Given a numpy array of shape (n, 2) representing the shapes of n images
#   (stored as numpy arrays), and 'trans', the shape (n, 2, 3) numpy array
#   used to represent affine transformations of n images in this script,
#   return a shape (n, 2) numpy array representing the shapes of n images
#   after applying .rotate(expand = True) to the PIL Image version of each
#   image
def rotate_shapes(shapes, trans):
    #   For each image shape, form a list of coordinates for the corners of the
    #   rectangle representing the image
    start_points = np.array(
        [
            [(0, 0), (0, shapes[i, 1]), (shapes[i, 0], 0), shapes[i]]
            for i in range(shapes.shape[0])
        ],
        dtype = np.float64
    )
    #   Apply the appropriate rotation to each of those points for each image
    end_points = np.array(
        [
            (trans[i, :, :2] @ start_points[i, :, :].T).T
            for i in range(shapes.shape[0])
        ],
        dtype = np.float64
    )
    #   Return the image dimensions of the rotated image by forming a rectangle
    #   that fits all 4 corners of each image
    end_shape = np.ceil(
        np.max(end_points, axis = 1) - np.min(end_points, axis = 1)
    ).astype(int)
    return end_shape
    
#   sample_info: pd.DataFrame of sample information with rows representing
#       samples and containing a 'raw_image_path' column with paths to full-
#       resolution images
#   theta: np.array (shape (N,)) containing rotations of every sample in radians
#   trans_img: np.array (shape (N, 2)) giving translations of every sample
#       (after any rotations)
#   max0, max1: integers giving shape of merged image to create
#
#   Return the combined full-resolution image containing all samples
#   after performing transformations (an np.array). Each sample occupies one
#   channel of the image (so RGB is reduced to grayscale). Intended for use
#   in creating the Samui-viewable image
def merge_image_fullres(sample_info, theta, trans_img, max0, max1):
    #   Initialize the combined tiff. Determine the boundaries by computing the
    #   maximal coordinates in each dimension of each rotated and translated
    #   image
    combined_img = np.zeros((max0 + 1, max1 + 1,6), dtype = np.uint8)
    weights = np.zeros((max0 + 1, max1 + 1, 6), dtype = np.uint8)
    
    for array in range(sample_info.shape[0]):
         #img_pil = Image.open(sample_info['raw_image_path'].iloc[array])
         #Rotate about the top left corner of the image
         theta_deg = 180 * theta[array] / np.pi # '.rotate' uses degrees, not radians
         #img = np.array(img_pil.rotate(theta_deg, expand = True, fillcolor=BACKGROUND_COLOR), dtype = np.uint16)
         
         imgRaw = tifffile.imread(sample_info['raw_image_path'].iloc[array])
         img_reshape = np.transpose(imgRaw, (1, 2, 0))
         
         #img = np.zeros(img_reshape.shape, dtype='uint8')
         for i in range(img_reshape.shape[2]):   
             img_pil = Image.fromarray(img_reshape[:,:,i])
             img_rotated = np.array(img_pil.rotate(theta_deg, expand = True, fillcolor=BACKGROUND_COLOR), dtype = np.uint8)
             if (i == 0): img = np.zeros((img_rotated.shape[0],img_rotated.shape[1], 6), dtype='uint8')
             img[:,:,i] = img_rotated
             
         #   "Place this image" on the combined image, considering translations.
         #   Use separate channels for each image
         #   "Place this image" on the combined image, considering translations.
         #   Use separate channels for each image
         # combined_img[
         #                trans_img[array, 0]: trans_img[array, 0] + img.shape[0],
         #                trans_img[array, 1]: trans_img[array, 1] + img.shape[1],
         #                :
         #            ] += img#.reshape((img.shape[0], img.shape[1], 6))
         #   Count how many times a pixel is added to
         
         combined_img[
             trans_img[array, 0]: trans_img[array, 0] + img.shape[0],
             trans_img[array, 1]: trans_img[array, 1] + img.shape[1],
             :
         ] = np.maximum(combined_img[
             trans_img[array, 0]: trans_img[array, 0] + img.shape[0],
             trans_img[array, 1]: trans_img[array, 1] + img.shape[1],
             :
         ], img)
         
         
         weights[
                 trans_img[array, 0]: trans_img[array, 0] + img.shape[0],
                 trans_img[array, 1]: trans_img[array, 1] + img.shape[1],
                 :
             ] += 1
    #   Fill in empty pixels with background color
    combined_img[weights[:, :, 0] == 0, :] = BACKGROUND_COLOR
    #   Average the color across all images that overlap a given pixel
    weights[weights == 0] = 1
    combined_img = (combined_img / weights).astype(np.uint8)
    return combined_img
    
################################################################################
#   Define the transformation matrix
################################################################################

#   Check if the adjusted estimates exist (i.e. if we already ran this script to
#   display in Samui, annotated ROIs, then ran '03-adjust_transform.py' to
#   refine the initial transformations from ImageJ). If not, use the initial
#   estimate from ImageJ
#if file_suffix == 'adjusted':
#    #   Read in the adjusted transformations
#    estimates_df = pd.read_csv(estimate_path)
#    estimates_df = estimates_df.loc[estimates_df['adjusted']]
#else:
estimates_df = sample_info.rename(
    {
        'initial_transform_x': 'x',
        'initial_transform_y': 'y',
        'initial_transform_theta': 'theta'
    },
    axis = 1
)

#   x and y switch, so angles invert
estimates_df['theta'] *= -1
theta = np.array(estimates_df['theta'])

#   Array defining an affine transformation:
#   [x0', y0'] = trans[0] @ [x0, y0, 1]
trans = np.array(
    [
        [
            [
                np.cos(estimates_df['theta'].iloc[i]),
                -1 * np.sin(estimates_df['theta'].iloc[i]),
                estimates_df['x'].iloc[i]
            ],
            [
                np.sin(estimates_df['theta'].iloc[i]),
                np.cos(estimates_df['theta'].iloc[i]),
                estimates_df['y'].iloc[i]
            ]
        ]
        for i in range(estimates_df.shape[0])
    ],
    dtype = np.float64
)

#   Flip x and y to follow Samui conventions. Translations must be an integer
#   number of pixels
trans[:, :, 2] = np.flip(np.round(trans[:, :, 2]), axis = 1)

################################################################################
#   Determine image shapes and initialize the combined image
################################################################################

#   Loop through all samples and grab image dimensions so we can compute
#   the dimensions of the combined image. Doesn't load the images into memory
img_shapes = []
for i in range(sample_info.shape[0]):
    tif = tifffile.TiffFile(sample_info['raw_image_path'].iloc[i])
    img_shapes.append(tif.pages[0].shape)
    tif.close()

img_shapes = np.array(img_shapes)
rotated_shapes = rotate_shapes(img_shapes[:, :2], trans)

trans_img = np.array(
    [
        trans[i, :, 2] +
        adjust_trans(np.flip(img_shapes[i, :2]), theta[i])
        for i in range(sample_info.shape[0])
    ],
    dtype = np.float64
)
trans_img = np.round(trans_img).astype(int)

#   Equally translate images and spot coords such that the minimum pixel
#   occupied by any image is at (0, 0) (also potentially saving memory)
trans[:, :, 2] -= np.min(trans_img, axis = 0)
trans_img -= np.min(trans_img, axis = 0)

#   For now, just read in one JSON file and assume all samples are roughly on
#   the same spatial scale (spots have the same diameter in pixels)
#json_path = os.path.join(sample_info['spaceranger_dir'].iloc[0], 'scalefactors_json.json')
#with open(json_path, 'r') as f:
#    spaceranger_json = json.load(f)

#m_per_px = SPOT_DIAMETER_JSON_M / spaceranger_json['spot_diameter_fullres']

################################################################################
#   Read in and translate images and spot coordinates
################################################################################

max0, max1 = np.max(rotated_shapes + trans_img, axis = 0)

#   Write the full-resolution combined TIFF needed by Samui
combined_img = merge_image_fullres(sample_info, theta, trans_img, max0, max1)
with tifffile.TiffWriter(img_out_fullres, bigtiff = True) as tiff:
    tiff.write(combined_img)

#   Write the low-resolution combined PNG needed for the SpatialExperiment
#   object, as well as a combined spaceranger-compatible JSON

#combined_img, lowres_sf = merge_image_export(sample_info, theta, trans_img, max0, max1)
#Image.fromarray(combined_img).save(img_out_export_path)

    #   Overwrite the lowres scale factor and delete the highres one (which
    #   isn't defined for our case). Then write
#spaceranger_json['tissue_lowres_scalef'] = lowres_sf
#spaceranger_json.pop('tissue_hires_scalef')
#with open(json_out_path, 'w') as f:
#    json.dump(spaceranger_json, f)

#   Handle the spot coordinates
tissue_positions_list = []
for i in range(sample_info.shape[0]):
    #   Read in the tissue positions file to get spatial coordinates. Index by
    #   barcode + sample ID
    pattern = re.compile(r'^tissue_positions(_list|)\.csv$')
    this_dir = sample_info['spaceranger_dir'].iloc[i]
    tissue_path = [
            Path(os.path.join(this_dir, x)) for x in os.listdir(this_dir)
            if pattern.match(x)
        ][0]
    if '_list' in tissue_path.name: 
        tissue_positions = pd.read_csv(
            tissue_path,
            header = None,
            # Note the switch of x and y
            names = ["in_tissue", "row", "col", "y", "x"],
            index_col = 0
        )
    else:
        tissue_positions = pd.read_csv(
            tissue_path,
            skiprows = 1,
            # Note the switch of x and y
            names = ["in_tissue", "row", "col", "y", "x"],
            index_col = 0
        )
    tissue_positions.index = tissue_positions.index + '_' + sample_info.index[i]
    #   Apply affine transform of coordinates
    tissue_positions[['y', 'x']] = (
        trans[i] @ 
        np.array(tissue_positions.assign(ones = 1)[['y', 'x', 'ones']]).T
    ).T
    tissue_positions_list.append(tissue_positions)

#   Also export tissue positions with SpatialExperiment-friendly colnames

tissue_positions_r = pd.concat(tissue_positions_list).rename(
        {
            'row': 'array_row', 'col': 'array_col', 'y': 'pxl_row_in_fullres',
            'x': 'pxl_col_in_fullres'
        },
        axis = 1
    )
tissue_positions_r.index.name = 'key'
tissue_positions_r.to_csv(tissue_out_path, index = True)

