import numpy as np
import h5py
from deepposekit.models import load_model
from videoreader import VideoReader
import leap_utils as lu
import xarray_behave as xb
from leap_utils import preprocessing
# from leap_utils.utils import iswin, islinux, ismac
from pathlib import Path
import xarray as xr
import pandas as pd
import logging
import zarr
# import dqefopt
import os
# from tqdm import tqdm


def deepposekit(tracksfilename: str, savename:str, modelname: str):
    datename = os.path.split(os.path.split(tracksfilename)[0])[1]
    root = os.getcwd()

    model_filename = modelname  #f'../snakemake-workflows/deepposekit/models/{modelname}'
    logging.info(f'using model from {model_filename}')
    model = load_model(model_filename)

    logging.info(f'assembling data for {datename}')
    dataset = xb.assemble(datename, root, dat_path='dat', res_path='res')
    logging.info(dataset)

    logging.info(f'will save to {savename}')

    vr = VideoReader(dataset.attrs['video_filename'])
    logging.info(vr)
    frame_numbers, frame_indices = np.unique(dataset.nearest_frame, return_index=True)
    frame_indices = frame_indices
    CENTER = 1
    box_centers = dataset.body_positions[frame_indices, :, CENTER, :].values
    nb_frames = len(frame_numbers)
    logging.info(f'loading boxes from {nb_frames} frames.')

    skeleton = pd.read_csv(os.path.dirname(modelname) + '/skeleton_initialized.csv')
    body_parts_skeleton = skeleton['name'].tolist()
    logging.info('body parts from skeleton: {}.'.format(body_parts_skeleton))
    body_parts = ['head', 'neck', 'front_left_leg', 'middle_left_leg', 'back_left_leg', 'front_right_leg', 'middle_right_leg', 'back_right_leg', 'thorax', 'left_wing', 'right_wing', 'tail']
    logging.info('overriding with: {}.'.format(body_parts))

    nb_flies = len(dataset.flies)
    nb_parts = len(body_parts)
    frame_start = min(frame_numbers)
    frame_stop = max(frame_numbers)
    batch_size = 100
    box_size = (64, 64) if 'scaled' in modelname else (96, 96)

    batch_idx = list(range(0, nb_frames, batch_size))
    nb_batches = len(batch_idx) - 1
    # batch_idx.append(int(frame_stop - frame_start - 2))

    poses = np.full((frame_stop, nb_flies, nb_parts, 2), np.nan)
    # poses[:] = np.nan
    poses_confidence = np.full((frame_stop, nb_flies, nb_parts, 1), np.nan)
    box_centers_from_batches = np.full((frame_stop, nb_flies, 2), np.nan)
    # box_centers_from_batches[:] = np.nan
    # pbar = tqdm(total=frame_stop, unit='frames')
    print('batch_idx', len(batch_idx))
    print('nb_batches', nb_batches)
    print(vr)
    print('nb_frames (vr)', len(vr))
    print('frame_numbers', len(frame_numbers))
    print('frame_indices', len(frame_indices))
    for batch_num in range(nb_batches):
        logging.info(f"PROCESSING BATCH {batch_num:03d}/{nb_batches} (frames {frame_numbers[batch_idx[batch_num]]}-{frame_numbers[batch_idx[batch_num + 1]]}).")
        batch_frame_numbers = list(range(frame_numbers[batch_idx[batch_num]], frame_numbers[batch_idx[batch_num + 1]]))
        box_centers_batch = box_centers[batch_idx[batch_num]:batch_idx[batch_num + 1]]

        boxes, fly_id, fly_frame = preprocessing.export_boxes(frames=[frame[..., :1] for frame in vr[batch_frame_numbers]],
                                                              box_centers=box_centers_batch,
                                                              box_size=box_size,
                                                              box_angles=None,
                                                              background=np.array(dataset.body_positions.attrs['background'])[..., np.newaxis])

        fly_frame = fly_frame + batch_frame_numbers[0]

        if 'bgsub' in modelname:
            boxes = np.clip(np.diff(boxes.astype(np.float16), axis=-1), 0, 1000).astype(np.uint8)
        elif 'bg' in modelname:
            boxes = boxes
        else:
            boxes = boxes[..., :1]
        # logging.info(f'estimating poses:')
        predictions = model.predict(boxes, batch_size=100, verbose=0)
        positions, confidence, _ = np.split(predictions, (2, 3), -1)
        box_centers_batch_flat = box_centers_batch.reshape((-1, 2))
        for pred, conf, box_center, fly, frame in zip(positions, confidence, box_centers_batch_flat, fly_id, fly_frame):
            poses[frame, fly, ...] = pred
            poses_confidence[frame, fly, ...] = conf
            box_centers_from_batches[frame, fly, ...] = box_center
    #     pbar.update(fly_frame[-1])
    # pbar.close()
    # convert to xarray and save to zarr
    poses = poses[..., ::-1]  # all data are in y/x - except for deepposekit - swap here to have order of coords like everywhere else
    poses = xr.DataArray(poses,
                         dims=['frames', 'flies', 'poseparts', 'coords'],
                         coords={'poseparts': body_parts,
                                 'coords': ['y', 'x']})
    poses_confidence = xr.DataArray(poses_confidence,
                                    dims=['frames', 'flies', 'poseparts', 'confidence'],
                                    coords={'poseparts': body_parts})
    box_centers = xr.DataArray(box_centers_from_batches,
                               dims=['frames', 'flies', 'coords'],
                               coords={'coords': ['y', 'x']})

    ds = xr.Dataset({"poses": poses, "poses_confidence": poses_confidence, "box_centers": box_centers},
                    attrs={'model_filename': model_filename,
                           'box_size': box_size,
                           'modelname': modelname,
                           'bgsub': 'bgsub' in modelname,
                           'video_filename': dataset.attrs['video_filename'],
                           })

    logging.info(ds)
    logging.info(f'Saving to {savename}.')
    with zarr.ZipStore(savename, mode='w') as zarr_store:
        ds.to_zarr(store=zarr_store, compute=True)


logging.basicConfig(level=logging.DEBUG)

from snakemake.shell import shell

params = snakemake.params[0][snakemake.rule]
extra = ''.join([f'--{k} {v} ' for k,v in params.items()])
os.umask(2) # is hopefully g+rwx - is ug=rwx,o=rx
for out in snakemake.output:
    deepposekit(snakemake.input.tracks, out, modelname=params['modelname'])
