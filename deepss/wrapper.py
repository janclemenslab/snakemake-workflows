"""Annotate vibrations with a network."""
import logging
import numpy as np
import defopt
import os
import dss.utils, dss.data, dss.models, dss.event_utils, dss.predict
import scipy.signal as ss
import h5py
from glob import glob
from typing import List

import flammkuchen


# move to cli module
def deepss(data_name: str, save_name: str, model_save_name: str, data_key: str = 'samples',
        event_class_name: str = 'vibration', nb_channels: int = 16, event_thres: float = 0.5, event_tol: float = 0.01):
    """[summary]

    Args:
        data_name (str): [description]
        save_name (str): [description]. Defaults to None.
        model_save_name (str): [description]
        data_key (str): [description]. Defaults to 'samples'.
        nb_channels (int): Number of channels to take from data file. Defaults to 16.
        event_class_name (str): [description]. Defaults to 'pulse'.
        event_thres (float): [description]. Defaults to 0.75.
        event_tol (float): [description]. Defaults to 0.01 seconds (10ms).

    Raises:
        ValueError: if data_name or save_name are of unknown type (allowed: wav, h5, zarr, npy/npz)
    """
    # load model
    logging.info(f'loading parameters for {model_save_name}')
    params = dss.utils.load_params(model_save_name)

    logging.info(f'loading data for {data_name}')
    with h5py.File(data_name, 'r') as f:
        x = f[data_key][..., :nb_channels]
        samplerate = f.attrs['rate']
        # channel_names = f.attrs['analog_chans_in']

    logging.info(f'   filtering')
    sos_bp = ss.butter(5, [50, 1000], 'bandpass', output='sos', fs=samplerate)
    x = ss.sosfiltfilt(sos_bp, x, axis=0).astype(np.float16)

    try:
        event_index = params['class_names'].index(event_class_name)
    except ValueError:
        logging.info(f'model does not predict events of the type "{event_class_name}".')
        event_index = None

    logging.info(f'   annotating {x.shape[0]/samplerate:1.2f} seconds')
    events, segments, class_probabilities = dss.predict.predict(x, model_save_name, params)
    logging.info(f'   saving')

    # save as xb
    logging.info(f'   saving results to "{save_name}".')

    # TODO make predict work with multiple event types (list of event_index) - not just vibrations
    os.makedirs(os.path.dirname(save_name), exist_ok=True)

    event_names = [k for k in events.keys() if k != 'samplerate_Hz']
    segment_names = [k for k in segments.keys() if k != 'samplerate_Hz' and k != 'noise']
    d = {'segment_names': segment_names,
        'segment_probabilities':  [segments[segment_name]['probabilities'] for segment_name in segment_names],
        'segment_labels': [segments[segment_name]['samples'] for segment_name in segment_names],
        'event_names': event_names,
        'event_probabilities': [events[event_name]['probabilities'] for event_name in event_names ],
        'event_indices': [events[event_name]['seconds'] * events['samplerate_Hz'] for event_name in event_names ],
        'samplerate_Hz': params['samplerate_y_Hz'],
        }
    flammkuchen.save(save_name, d)


from snakemake.shell import shell

logging.basicConfig(level=logging.INFO)
params = snakemake.params[0][snakemake.rule]
for out in snakemake.output:
    deepss(snakemake.input[0], save_name=out, model_save_name=params['modelname'])
