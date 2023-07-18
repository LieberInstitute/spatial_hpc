import os
import numpy as np
from cellpose import models, io
from cellpose.io import imread
import pyhere
import pandas as pd
from pathlib import Path

#   Pretrained model based on 'cyto' and iteratively refined in the GUI
model_path = pyhere.here(
    'processed-data', 'spot_deconvo', 'cellpose', 'richard', 'models',
    'CP_20220415_152031'
)
cell_diameter = None
channel = 2 # DAPI

#   Path to excel sheet containing sample info; directory containing .tif
#   Visium-IF images
sample_info_path = pyhere.here('raw-data', 'sample_info', 'VSPG_HPC_20230626.xlsx')
img_dir = pyhere.here('processed-data', 'Images', 'VistoSeg', 'VSPG')
mask_dir = pyhere.here('processed-data', 'spot_deconvo', 'cellpose', 'masks')

Path(mask_dir).mkdir(parents=True, exist_ok=True)

#   Determine paths to IF images; read in just the one for this sample
sample_info = pd.read_excel(sample_info_path)[1:9]
sample_info = sample_info.iloc[[2,7]]
sample_ids = sample_info['Slide #'] + '_' + sample_info['Array #']
sample_id = sample_ids[int(os.environ['SGE_TASK_ID']) - 1]

print(f'Segmenting DAPI for sample {sample_id}.')

img = imread(pyhere.here(img_dir, sample_id + '.tif'))[channel, :, :]

#   Load the pretrained model and process images using it
model = models.CellposeModel(
    gpu = True, pretrained_model = model_path.as_posix()
)
masks, flows, styles = model.eval(img, diameter=cell_diameter)

#   Save PNG version of the masks to visually inspect results
mask_png = str(pyhere.here(mask_dir, f'{sample_id}_mask.png'))
io.save_to_png(img, masks, flows, mask_png)

#   Save masks
mask_npy = str(
    pyhere.here(mask_dir, f'{sample_id}_mask.npy')
)
np.save(mask_npy, masks)