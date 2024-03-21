from pyhere import here
from pathlib import Path
import session_info
import os
import re
import json

import pandas as pd
import numpy as np

out_path = here('processed-data', '10-image_stitching', 'sample_info_clean.csv')

sample_info_path = Path(
    here('raw-data', 'sample_info', 'samSheet.xlsx')
)

xml_map_path = here('raw-data', 'sample_info', 'sample_key.csv')


# for this project this should be directly inside 01_spaceranger directory
#sRanger = "spaceranger_2022-04-12_SPag033122"
#sRanger = "spaceranger_2023-01-31_round9"
# symlink directory for spaceranger output for samui donors
sRanger = "spaceranger-all"

################################################################################
#   Functions
################################################################################

#   Given the contents of the spaceranger '_invocation' file as a single string
#   'file_string', return the path to the image used as input to spaceranger
def parse_sr_invocation(file_string):
    #   Find the line containing the image path, and confirm there is only
    #   one image 
    match = re.findall(r'tissue_image_paths.*=.*\[".*"\]', file_string)
    assert len(match) == 1

    #   Grab just the path and return
    path = Path(re.sub(r'(tissue_image_paths|[ "\[\]=,])', '', match[0]))
    assert path.exists()
    return path

#   Get a numpy array representing transformation matrices from ImageJ (of shape
#   (n, 2, 3) for n capture areas), return theta (of shape (n,)), the angles
#   determined by the rotation-matrix portion of 'mat'
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
sample_info = (
    pd.read_excel(sample_info_path).loc[:, ['Brain', 'Slide #', 'Array #']]
)

#   Make brain number consistent
sample_info['Brain'] = (
    sample_info['Brain']
        .astype(str)
        #   Remove additional pieces like '-left' or 'Hs_' from brain number
        .replace(to_replace = r'^.*(Br[0-9]{4}).*', value = r'\1', regex = True)
)

#   Make slide number consistent
sample_info['Slide #'] = (
    sample_info['Slide #']
        #   One slide has a mistaken '-' that must be removed
        .replace(r'V([0-9]{2})-([A-Z][0-9]{2})', r'V\1\2', regex = True)
)

sample_info.index = sample_info['Slide #'] + '_' + sample_info['Array #']

################################################################################
#   Merge in info about where initial transforms are
################################################################################

xml_map = pd.read_csv(xml_map_path, index_col = 'Slide')
xml_map.index.name = 'sample_id'
xml_map = xml_map.loc[:, ['Brain', 'XML file name']].copy()
xml_map['Slide #'] = [x.split('_')[0] for x in xml_map.index]
xml_map['Array #'] = [x.split('_')[1] for x in xml_map.index]
#xml_map['In analysis'] = xml_map['In analysis'] == 'Yes'

#   Merge sample_info and xml_map by index, keeping the union of their columns.
#   There should be a more elegant way to do this...
sample_info = pd.merge(
    sample_info, xml_map, how = "outer", left_index = True, right_index = True
)
sample_info['Brain'] = sample_info['Brain_x'].combine_first(sample_info['Brain_y'])
sample_info['Slide #'] = sample_info['Slide #_x'].combine_first(sample_info['Slide #_y'])
sample_info['Array #'] = sample_info['Array #_x'].combine_first(sample_info['Array #_y'])
sample_info.drop(
    ['Brain_x', 'Slide #_x', 'Brain_y', 'Slide #_y', 'Array #_x', 'Array #_y'],
    axis = 1, inplace = True
)

################################################################################
#   Determine the path to the full-resolution/raw images and spaceranger output
#   directories
################################################################################

sample_info['spaceranger_dir'] = [
    Path(x)
    for x in here(
        'processed-data', '01_spaceranger', sRanger, sample_info.index, 'outs',
        'spatial'
    )
]

#   Grab the full-resolution ("raw") image paths by extracting which images were
#   used as input to spaceranger
raw_image_paths = []
for spaceranger_dir in sample_info['spaceranger_dir']:
    invocation_path = spaceranger_dir.parent.parent / '_invocation'
    
    #   Spaceranger has not been run for all samples
    if invocation_path.exists():
        with open(invocation_path, 'r') as f:
            invocation_str = f.read()
        
        raw_image_paths.append(parse_sr_invocation(invocation_str))
    else:
        raw_image_paths.append(None)

sample_info['raw_image_path'] = raw_image_paths

################################################################################
#   Add the initial transformation estimates for image stitching from ImageJ
################################################################################

transform_df_list = []
for imagej_xml_path in sample_info['XML file name'].dropna().unique():
    #   Open the ImageJ XML output
    with open(here(imagej_xml_path)) as f:
        imagej_xml = f.read()

    #   Clean the file; make sure new lines only separate XML elements
    imagej_xml = re.sub('\n', '', imagej_xml)
    imagej_xml = re.sub('\>', '>\n', imagej_xml)

    array_nums = [
        re.sub('[_/]', '', x)
        for x in re.findall(r'_[ABCD]1/', imagej_xml)
    ]
    slide_nums = re.findall(r'V[0-9]{2}[A-Z][0-9]{2}-[0-9]{3}', imagej_xml)

    #   Some older samples use donor + array instead of slide + array as sample IDs.
    #   It so happens that those samples are always on one slide, so this ugly
    #   workaround grabs the correct slide for each such sample while also correctly
    #   handling new samples that use slide + array. It takes advantage of the fact
    #   that for old samples, the slide number is always in the title of the XML
    #   project, which the below regexs look for
    #if len(slide_nums) == len(re.findall(r'title="V[0-9]{2}[A-Z][0-9]{2}', imagej_xml)):
    #    #   Just use the slide in the title
    #    slide_nums = slide_nums * len(array_nums)
    #elif len(re.findall(r'title="V[0-9]{2}[A-Z][0-9]{2}', imagej_xml)) == 1:
    #    #   Otherwise use the slides in the image paths, not the title
    #slide_nums = slide_nums[1:]
    
    #   Grab the transformation matrices and import as numpy array with the
    #   structure used in 01-samui_test.py
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

    #   Grab sample info for this slide, ordered how the array numbers
    #   appear in the ImageJ output
    this_sample_info = sample_info.loc[
        [slide_nums[i] + '_' + array_nums[i] for i in range(len(array_nums))]
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
                'Array #': array_nums,
                'Slide #': slide_nums
            }
        )
    )

#   Combine across donors, then merge into sample_info
transform_df = pd.concat(transform_df_list)
index = sample_info.index
sample_info = pd.merge(
    sample_info, transform_df, how = 'left', on = ['Slide #', 'Array #']
)
sample_info.index = index

sample_info.to_csv(out_path)

session_info.show()
