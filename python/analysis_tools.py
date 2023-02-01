
# Import MNE-NIRS processing
import mne
from mne_nirs.channels import get_long_channels
from mne_nirs.channels import picks_pair_to_idx
from mne_nirs.datasets import fnirs_motor_group
from mne.preprocessing.nirs import beer_lambert_law
from mne_nirs.signal_enhancement import enhance_negative_correlation

# Import MNE processing
from mne.viz import plot_compare_evokeds
from mne import Epochs, events_from_annotations, set_log_level
from mne.filter import resample

# Other Tooling
import pandas as pd
import numpy as np
import os
from collections import defaultdict
from copy import deepcopy
from itertools import compress
import matplotlib.pyplot as plt
import importlib


# Local functions
import quality_eval

importlib.reload(quality_eval)

def individual_analysis(bids_path, trigger_id, variable_epoch_time):
    """
    TLDR:
        This function takes in the file path to the BIDS directory and the dictionary that renames numeric triggers.
        This function returns raw haemodynamic (from the beer lambert function) data per MNE specs and epochs for the triggers.

    Examples of each variables:
        bids_path (str) = '../../LabResearch/IndependentStudy/DataAnalysis'
        trigger_id (dict) = {'4': 'Control', '2': 'Neutral', '3': 'Inflammatory', '1':'Practice'}

    Documentation:
        raw_haemo: https://mne.tools/stable/auto_tutorials/preprocessing/70_fnirs_processing.html#sphx-glr-auto-tutorials-preprocessing-70-fnirs-processing-py
        epochs: https://mne.tools/stable/generated/mne.Epochs.html#mne-epochs

    """
    # Read data with annotations in BIDS format
    # raw_intensity = read_raw_bids(bids_path=bids_path, verbose=False)
    raw_intensity = mne.io.read_raw_snirf(bids_path, verbose=True, preload=False)

    raw_intensity = get_long_channels(raw_intensity, min_dist=0.01)
    
    # Rename the numeric triggers for ease of processing later
    raw_intensity.annotations.rename(trigger_id)

    raw_od = quality_eval.signal_preprocessing(raw_intensity)

    # Convert to haemoglobin and filter
    raw_haemo = beer_lambert_law(raw_od, ppf=0.1)
    raw_haemo = raw_haemo.filter(0.02, 0.3,
                                 h_trans_bandwidth=0.1, l_trans_bandwidth=0.01,
                                 verbose=False)

    # Apply further data cleaning techniques and extract epochs
    raw_haemo = enhance_negative_correlation(raw_haemo)
    # Extract events but ignore those with
    # the word Ends (i.e. drop ExperimentEnds events)

    # TODO: Need to add something here to be able to work in custom triggers
    events, event_dict = events_from_annotations(raw_haemo, verbose=False)

    if variable_epoch_time:
        epochs = dynamic_time_epoch_generation(raw_haemo, event_dict, events)
    else:
        # Remove all STOP triggers to hardcode duration to 30 secs per MNE specs
        #TODO: We'll need to remove this for all other datasets
        events = events[::2]
        
        epochs = Epochs(raw_haemo, events, event_id=event_dict, tmin=-1, tmax=15,
                        reject=dict(hbo=200e-6), reject_by_annotation=True,
                        proj=True, baseline=(None, 0), detrend=0,
                        preload=True, verbose=False)
        # Doing this to ensure that if variable_epoch_time is True or False the format
        # of the data remains the same.
        epochs = [epochs]

    print(epochs)

    return epochs


def aggregate_epochs(paths, trigger_id, variable_epoch_time):
    """
    TLDR:
        Cycles through the participants in bids folders and returns epochs based on the trigger associated with it
        This function takes in the root directory, trigger_id dict and ignore parameters.
        This function returns a dict with the keys being trigger names and the value being an array of epochs.

    Examples of each variables:
        root_dir ([str]) = ['../../LabResearch/IndependentStudy/DataAnalysis'] generated from import_data_folder function
        trigger_id (dict) = {'4': 'Control', '2': 'Neutral', '3': 'Inflammatory', '1':'Practice'}

    Documentation:
        epochs: https://mne.tools/stable/generated/mne.Epochs.html#mne-epochs

    """
    all_epochs = defaultdict(list)

    for f_path in paths:

        epochs = individual_analysis(f_path, trigger_id, variable_epoch_time)
        for epoch in epochs:
            for cidx, condition in enumerate(epoch.event_id):
                all_epochs[condition].append(epoch[condition])
    
    return all_epochs


def extract_all_amplitudes(all_epochs, tmin, tmax):
    """
        TLDR:
            Takes in all_epochs dict and returns a data frame of all the the optical measurements taken during the experiment.
            The dataframe columns represent timestamps and the rows represent hbo/hbr of each users fNIRs session.

        Examples of each variables:
            all_epochs (dict) = 
                    ({'Control': [<Epochs |  3 events (all good), -1.25 - 15 sec, baseline -1.25 – 0 sec, ~171 kB, data loaded,
                    'Control': 3>,
                    <Epochs |  3 events (all good), -1.25 - 15 sec, baseline -1.25 – 0 sec, ~171 kB, data loaded,
                    'Control': 3>,
                    <Epochs |  3 events (all good), -1.25 - 15 sec, baseline -1.25 – 0 sec, ~171 kB, data loaded,
                    'Control': 3>])
            columns (array) = ['Name 1', 'Name 2', 'Name 3'] 
                * note: array length must match the number of columns
            tmin (int/float) = In relation to 0 (time of trigger) what is the minimum time from that point you want to analyze (can be negative)
            tmax (inf/float) = In relation to 0 (time of trigger) what is the maximum time from that point you want to analyze

        Documentation:
            epochs: https://mne.tools/stable/generated/mne.Epochs.html#mne-epochs

    """
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


def extract_average_amplitudes(all_epochs, tmin, tmax):
    """
        TLDR:
            Takes in all_epochs dict and returns a data frame of the average measurements taken during the experiment per subject and condition.
            The dataframe columns represent timestamps and the rows represent hbo/hbr of each users fNIRs session.

        Examples of each variables:
            all_epochs (dict) = 
                    ({'Control': [<Epochs |  3 events (all good), -1.25 - 15 sec, baseline -1.25 – 0 sec, ~171 kB, data loaded,
                    'Control': 3>,
                    <Epochs |  3 events (all good), -1.25 - 15 sec, baseline -1.25 – 0 sec, ~171 kB, data loaded,
                    'Control': 3>,
                    <Epochs |  3 events (all good), -1.25 - 15 sec, baseline -1.25 – 0 sec, ~171 kB, data loaded,
                    'Control': 3>])
            tmin (int/float) = In relation to 0 (time of trigger) what is the minimum time from that point you want to analyze (can be negative)
            tmax (inf/float) = In relation to 0 (time of trigger) what is the maximum time from that point you want to analyze

        Documentation:
            epochs: https://mne.tools/stable/generated/mne.Epochs.html#mne-epochs

    """
    columns = ['ID', 'Chroma', 'Condition', 'Value']
    df = pd.DataFrame(columns=columns)
    temporal_measurements = []
    columns = ['ID', 'Chroma', 'Condition', 'Value']

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


def dynamic_time_epoch_generation(raw_haemo, event_dict, events):
    # print(len(events))
    # print(raw_haemo)

    for index, event in enumerate(events):
        if index % 2 == 0 and index != 0:
            prev_event_time = events[index - 2][0]
            current_event_time = events[index - 1][0]
            task_len = current_event_time - prev_event_time

            epochs = Epochs(raw_haemo, events, event_id=event_dict, tmin=-1, tmax=15,
                        reject=dict(hbo=200e-6), reject_by_annotation=True,
                        proj=True, baseline=(None, 0), detrend=0,
                        preload=True, verbose=False)        

    return epochs