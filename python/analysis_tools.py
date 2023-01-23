
# Import MNE-NIRS processing
import mne
from mne_nirs.channels import get_long_channels
from mne_nirs.channels import picks_pair_to_idx
from mne_nirs.datasets import fnirs_motor_group
from mne.preprocessing.nirs import beer_lambert_law, optical_density, temporal_derivative_distribution_repair, scalp_coupling_index
from mne_nirs.signal_enhancement import enhance_negative_correlation

# Import MNE processing
from mne.viz import plot_compare_evokeds
from mne import Epochs, events_from_annotations, set_log_level

# Other Tooling
import pandas as pd
import numpy as np
import os
from collections import defaultdict
from copy import deepcopy
from itertools import compress


def individual_analysis(bids_path, trigger_id):

    # Read data with annotations in BIDS format
    # raw_intensity = read_raw_bids(bids_path=bids_path, verbose=False)
    raw_intensity = mne.io.read_raw_snirf(bids_path, verbose=True, preload=False)
    raw_intensity = get_long_channels(raw_intensity, min_dist=0.01)
    
    channel_types = raw_intensity.copy()
    # print(channel_types)
    
    raw_intensity.annotations.rename(trigger_id)

    # Convert signal to optical density and determine bad channels
    raw_od = optical_density(raw_intensity)
    sci = scalp_coupling_index(raw_od, h_freq=1.35, h_trans_bandwidth=0.1)
    raw_od.info["bads"] = list(compress(raw_od.ch_names, sci < 0.5))
    # raw_od.interpolate_bads()

    # Downsample and apply signal cleaning techniques
    raw_od.resample(0.8)
    raw_od = temporal_derivative_distribution_repair(raw_od)

    # Convert to haemoglobin and filter
    raw_haemo = beer_lambert_law(raw_od, ppf=0.1)
    raw_haemo = raw_haemo.filter(0.02, 0.3,
                                 h_trans_bandwidth=0.1, l_trans_bandwidth=0.01,
                                 verbose=False)

    # Apply further data cleaning techniques and extract epochs
    raw_haemo = enhance_negative_correlation(raw_haemo)
    # Extract events but ignore those with
    # the word Ends (i.e. drop ExperimentEnds events)
    events, event_dict = events_from_annotations(raw_haemo, verbose=False)
    
    # Remove all STOP triggers to hardcode duration to 30 secs per MNE specs
    events = events[::2]
    # print(events)

    epochs = Epochs(raw_haemo, events, event_id=event_dict, tmin=-1, tmax=15,
                    reject=dict(hbo=200e-6), reject_by_annotation=True,
                    proj=True, baseline=(None, 0), detrend=0,
                    preload=True, verbose=False)

    return raw_haemo, epochs


def aggregate_epochs(root_dir, trigger_id, ignore):
    all_evokeds = defaultdict(list)

    subjects = os.listdir(f'{root_dir}/BIDS_Anon/')

    for sub in subjects:
        if sub not in ignore:
            # Create path to file based on experiment info
            f_path = f'{root_dir}/BIDS_Anon/{sub}/nirs/{sub}_task-AnonCom_nirs.snirf'

            # Analyze data and return both ROI and channel results
            raw_haemo, epochs = individual_analysis(f_path, trigger_id)

            for cidx, condition in enumerate(epochs.event_id):
                # all_evokeds[condition].append(epochs[condition].average())
                all_evokeds[condition].append(epochs[condition])

    return all_evokeds


def extract_all_amplitudes(all_epochs, columns, tmin, tmax):
    df = pd.DataFrame(columns=columns)
    temporal_measurements = []

    for idx, epoch in enumerate(all_epochs):
        subj_id = 0
        for subj_data in all_epochs[epoch]:
            subj_id += 1
            # can be either "hbo", "hbr", or both
            for chroma in ["hbo", "hbr"]:
                data = deepcopy(subj_data.average(picks=chroma))
                value = data.crop(tmin=tmin, tmax=tmax).data * 1.0e6
                
                # Reshape the data to be a flat numpy array
                value = np.reshape(value, -1)
                temporal_measurements.append(value)

    temporal_measurements = np.array(temporal_measurements)
    measurement_df = pd.DataFrame(temporal_measurements)

    return measurement_df


def extract_average_amplitudes(all_epochs, columns, tmin, tmax):
    df = pd.DataFrame(columns=columns)
    temporal_measurements = []

    for idx, epoch in enumerate(all_epochs):
        subj_id = 0
        for subj_data in all_epochs[epoch]:
            subj_id += 1
            # can be either "hbo", "hbr", or both
            for chroma in ["hbo", "hbr"]:
                data = deepcopy(subj_data.average(picks=chroma))
                value = data.crop(tmin=tmin, tmax=tmax).data * 1.0e6
                
                # Reshape the data to be a flat numpy array
                value = np.reshape(value, -1)
                temporal_measurements.append(value)

                # Placeholder while we see if PCA gives better results (also how the MNE tutorial said to do it though)
                avg_val = data.crop(tmin=tmin, tmax=tmax).data.mean() * 1.0e6

                # Append metadata and extracted feature to dataframe
                this_df = pd.DataFrame(
                    {'ID': subj_id, 'Chroma': chroma, 'Condition': epoch, 'Value': avg_val}, index=[0])
                df = pd.concat([df, this_df], ignore_index=True)

    df['Value'] = pd.to_numeric(df['Value'])  # some Pandas have this as object
    return df