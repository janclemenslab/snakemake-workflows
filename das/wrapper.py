"""Annotate vibrations with a network."""
import logging
import numpy as np
import os
import das.utils, das.event_utils, das.predict, das.annot
import scipy.signal as ss
from typing import List, Union
import xarray_behave as xb


def sampletimes_to_timestamps(events, segments, ds, start_index: int = 0, suffix: str = ''):
    event_times = das.annot.Events()
    fs_song = ds.song_raw.attrs['sampling_rate_Hz']

    # EVENTS
    detected_event_names = []
    if 'sequence' in events:
        detected_event_names = np.unique(events['sequence'])

    if len(detected_event_names) > 0 and detected_event_names[0] is not None:
        logging.info(f"   adding {detected_event_names} to annotations.")
        for name in detected_event_names:
            event_times.add_name(name, category='event')

        logging.info(f"   found {len(events['seconds'])} instances of events '{detected_event_names}'.")
        event_samples = (np.array(events['seconds']) * fs_song + start_index).astype(np.uintp)

        # make sure all detected events are within bounds
        event_samples = event_samples[event_samples >= 0]
        event_samples = event_samples[event_samples < ds.sampletime.shape[0]]

        event_seconds = ds.sampletime[event_samples]
        for name_or_index, seconds in zip(events['sequence'], event_seconds):
            if type(name_or_index) is int:
                event_name = str(events['names'][name_or_index])
            else:
                event_name = str(name_or_index)
            event_times.add_time(event_name + str(suffix), seconds, seconds)  #, category='event')

    # SEGMENTS
    detected_segment_names = []
    if 'sequence' in segments:
        detected_segment_names = np.unique(segments['sequence'])
        # if these are indices, get corresponding names
        if len(detected_segment_names) and type(detected_segment_names[0]) is not str and type(
                detected_segment_names[0]) is not np.str_:
            detected_segment_names = [segments['names'][ii] for ii in detected_segment_names]

    if len(detected_segment_names) > 0:  # and detected_segment_names[0] is not None:
        logging.info(f"   adding {detected_segment_names} to annotations.")
        for name in detected_segment_names:
            event_times.add_name(name, category='segment')

        logging.info(f"   found {len(segments['onsets_seconds'])} instances of segments '{detected_segment_names}'.")
        onsets_samples = (np.array(segments['onsets_seconds']) * fs_song + start_index).astype(np.uintp)
        offsets_samples = (np.array(segments['offsets_seconds']) * fs_song + start_index).astype(np.uintp)

        # make sure all detected segments are within bounds
        onsets_samples = onsets_samples[onsets_samples >= 0]
        offsets_samples = offsets_samples[offsets_samples >= 0]
        onsets_samples = onsets_samples[onsets_samples < ds.sampletime.shape[0]]
        offsets_samples = offsets_samples[offsets_samples < ds.sampletime.shape[0]]

        onsets_seconds = ds.sampletime[onsets_samples]
        offsets_seconds = ds.sampletime[offsets_samples]
        for name_or_index, onset_seconds, offset_seconds in zip(segments['sequence'], onsets_seconds, offsets_seconds):
            if type(name_or_index) is not str and type(detected_segment_names[0]) is not np.str_:
                segment_name = segments['names'][name_or_index]
            else:
                segment_name = str(name_or_index)

            event_times.add_time(segment_name + str(suffix), onset_seconds, offset_seconds)  #, category='segment')

    return event_times


def run_das(data_name: str, save_name: str, model_save_name: Union[List[str], str], as_timestamps: bool = True):
    """[summary]

    Args:
        data_name (str): [description]
        save_name (str): [description]. Defaults to None.
        model_save_name (str): [description]
        data_key (str): [description]. Defaults to 'samples'.
        nb_channels (int): Number of channels to take from data file. Defaults to 16.
        event_thres (float): [description]. Defaults to 0.5.
        event_tol (float): [description]. Defaults to 0.01 seconds (10ms).

    Raises:
        ValueError: if data_name or save_name are of unknown type (allowed: wav, h5, zarr, npy/npz)
    """

    logging.info(f'loading data for {data_name}')
    date_name = os.path.basename(os.path.dirname(data_name))
    ds = xb.assemble(date_name, include_poses=False, include_tracks=False, resample_video_data=False)

    samplerate = ds.song_raw.attrs['sampling_rate_Hz']
    logging.info(f'   filtering')
    sos_bp = ss.butter(5, [50, 1000], 'bandpass', output='sos', fs=samplerate)
    x = ss.sosfiltfilt(sos_bp, ds.song_raw, axis=0).astype(np.float16)

    if type(model_save_name) == str:
        model_save_name = [model_save_name]

    events = {'seconds': np.empty((0,)), 'sequence': np.empty((0,)), 'names': np.empty((0,))}
    segments = {
        'onsets_seconds': np.empty((0,)),
        'offsets_seconds': np.empty((0,)),
        'sequence': np.empty((0,)),
        'names': np.empty((0,)),
    }

    for model_name in model_save_name:
        # load model
        logging.info(f'loading parameters for {model_name}')
        model, params = das.utils.load_model_and_params(model_name)

        logging.info(f'   annotating {x.shape[0]/samplerate:1.2f} seconds')
        this_events, this_segments, _, _ = das.predict.predict(x, model=model, params=params, verbose=2)

        if this_events:
            for key in events.keys():
                events[key] = np.append(events[key], this_events[key])
        if this_segments:
            for key in segments.keys():
                if key == 'names':  # remove the noise (=no song) segment
                    this_segments[key] = [k for k in this_segments[key] if k != 'noise']
                segments[key] = np.append(segments[key], this_segments[key])

    # as timestamps
    evt = sampletimes_to_timestamps(events, segments, ds).to_df()
    # as samples
    evt_samples = das.annot.Events.from_predict(events, segments).to_df()
    evt['start_samples'] = (evt_samples['start_seconds'] * samplerate).astype(int)
    evt['stop_samples'] = (evt_samples['stop_seconds'] * samplerate).astype(int)

    logging.info(f"   Saving results to {save_name}.")
    os.makedirs(os.path.dirname(save_name), exist_ok=True)
    evt.to_csv(save_name)
    logging.info("Done.")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    print("snakemake.rule", snakemake.rule)
    print("snakemake.params[0]", snakemake.params[0])
    print("len(snakemake.params)", len(snakemake.params))
    print([p for p in snakemake.params])
    params = snakemake.params[0][snakemake.rule]
    print(params)
    for out in snakemake.output:
        print(snakemake.input[0], out, params['modelname'])
        run_das(snakemake.input[0], save_name=out, model_save_name=params['modelname'])

    # # # LOCAL TEST
    # dn = 'localhost-20221209_135308'
    # d = f'dat/{dn}/{dn}_daq.h5'
    # m = ['train_sep.res/pulse_20221222_101109', 'train_sep.res/sine_20221222_101035']
    # run_das(d, save_name=f'./{dn}_test.csv', model_save_name=m)
