import os
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')

from pyhere import here
from pathlib import Path
import session_info
import re
import json

import pandas as pd
import numpy as np

out_path = here('processed-data', 'VSPG_image_stitching', 'transformations.csv')
sample_info_path = Path(here('raw-data', 'sample_info', 'VSPG_HPC_20230626.xlsx'))

################################################################################
#   Functions
################################################################################

def theta_from_mat(mat):
    #   Determine theta by using arccos. Here the real angle might by theta_cos
    #   or -1 * theta_cos
    theta_cos = np.arccos(mat[:, 0, 0])
    #   Use 'theta_cos', but negate it when sin(theta) is negative, to account
    #   for the convention that arccos returns positive angles
    return theta_cos * (2 * (mat[:, 1, 0] > 0) - 1)


################################################################################
#   Clean sample info
################################################################################

#   Read in the sample sheet
sample_info = (pd.read_excel(sample_info_path).loc[1:9, ['Brain', 'Slide #', 'Array #']])
sample_info.index = sample_info['Slide #'] + '_' + sample_info['Array #']
#sample_info['xml_path'] = sample_info['Slide #'].apply(lambda x: 'processed-data/VSPG_image_stitching/335.xml' if x.endswith('335') else '')
sample_info['xml_path'] = sample_info['Slide #'].apply(lambda x: 'processed-data/VSPG_image_stitching/V12D07-335.xml' if x.endswith('335') else 'processed-data/VSPG_image_stitching/V12D07-332.xml')

################################################################################
#   Determine the path to the full-resolution/raw images and spaceranger output directories
################################################################################

sample_info['spaceranger_dir'] = [
    Path(x)
    for x in here(
        'processed-data', '01_spaceranger', 'spaceranger_VSPG', sample_info.index+'_'+sample_info.Brain, 'outs', 'spatial'
    )
]

sample_info['raw_image_path'] = [
    Path(x)
    for x in here(
        'processed-data', 'Images', 'VistoSeg', 'VSPG', sample_info.index+'.tif'
    )
]

################################################################################
#   Add the initial transformation estimates for image stitching from ImageJ
################################################################################

transform_df_list = []
for imagej_xml_path in sample_info['xml_path'].dropna().unique():
    #   Open the ImageJ XML output
    with open(here(imagej_xml_path)) as f:
        imagej_xml = f.read()
    arrays = re.findall(r'[ABCD]1', imagej_xml)
    unique_set = set()
    array_nums = [x for x in arrays if x not in unique_set and not unique_set.add(x)]
    
    slide_nums = re.findall(r'V[0-9]{2}[A-Z][0-9]{2}-[0-9]{3}', imagej_xml)
    #slide_nums = list(set(re.findall(r'33[0-9]_[ABCD]1', imagej_xml)))
    #   Clean the file; make sure new lines only separate XML elements
    imagej_xml = re.sub('\n', '', imagej_xml)
    imagej_xml = re.sub('\>', '>\n', imagej_xml)
    #   Grab the transformation matrices and import as numpy array with the structure used in 02-create_samui.py
    matrices = re.findall(r'matrix\(.*\)".*file_path=', imagej_xml)
    trans_mat = [re.sub(r'matrix\((.*)\).*', '\\1', x).split(',') for x in matrices]
    trans_mat = np.transpose(
        np.array(
            [[float(y) for y in x] for x in trans_mat], dtype = np.float64
        )
            .reshape((-1, 3, 2)),
        axes = [0, 2, 1]
    )
    assert trans_mat.shape[1:3] == (2, 3), "Improperly read ImageJ XML outputs"
    #   Grab sample info for this slide, ordered how the array numbers appear in the ImageJ output
    this_sample_info = sample_info.loc[
            [slide_nums[0] + '_' + array_nums[i] for i in range(len(array_nums))]
        ]
    #   Adjust translations to represent pixels in full resolution
    for i in range(len(array_nums)):
        #   Grab high-res scale factor and scale translations accordingly
        json_path = os.path.join(
            this_sample_info['spaceranger_dir'].iloc[i],
            'scalefactors_json.json'
        )
        with open(json_path, 'r') as f:
            spaceranger_json = json.load(f)
        trans_mat[i, :, 2] /= spaceranger_json['tissue_hires_scalef']
    #   Compile transformation information for this slide in a DataFrame
    transform_df_list.append(
        pd.DataFrame(
            {
                'initial_transform_x': trans_mat[:, 0, 2],
                'initial_transform_y': trans_mat[:, 1, 2],
                'initial_transform_theta': theta_from_mat(trans_mat),
                'Array #': this_sample_info['Array #'],
                'Slide #': this_sample_info['Slide #']
                
            }
        )
    )

#   Combine across donors, then merge into sample_info
transform_df = pd.concat(transform_df_list)
index = sample_info.index
sample_info = pd.merge(sample_info, transform_df, how = 'left', on = ['Slide #', 'Array #'])
sample_info.index = index

sample_info.to_csv(out_path)

session_info.show()
