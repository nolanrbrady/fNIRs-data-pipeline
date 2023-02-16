
# Import MNE-NIRS processing
import mne
from mne_nirs.channels import get_long_channels
from mne_nirs.channels import picks_pair_to_idx
from mne_nirs.datasets import fnirs_motor_group
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
import dynamic_interval_tools

importlib.reload(quality_eval)
importlib.reload(dynamic_interval_tools)

def individual_analysis(bids_path, trigger_id, variable_epoch_time, custom_triggers = False):
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

    raw_haemo = quality_eval.signal_preprocessing(raw_intensity)

    # Apply further data cleaning techniques and extract epochs
    raw_haemo = enhance_negative_correlation(raw_haemo)

    
    if custom_triggers:
        # TODO: Need to add something here to be able to work in custom triggers
        print('We need to add code to handle custom triggers')
    else:
        events, event_dict = events_from_annotations(raw_haemo, verbose=False)

    # Logic splits here since there are fundamental differences in how we handle Epoch
    # generation in dynamic intervals instead of block intervals.
    if variable_epoch_time:
        epochs = dynamic_interval_tools.epoch_generation(raw_haemo, event_dict, events)
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

    return epochs, raw_haemo, raw_intensity, bids_path


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
    # Temporary storage for the items we're using in all_data_df
    all_data = []

    columns = ['raw_haemo'] # 'epoch', 'condition', 

    for f_path in paths:

        epochs, raw_haemo, raw_intensity, path = individual_analysis(f_path, trigger_id, variable_epoch_time)
        for epoch in epochs:
            for cidx, condition in enumerate(epoch.event_id):
                all_epochs[condition].append(epoch[condition])
                epoch_data = {
                    'epoch': epoch,
                    'condition': condition,
                    'raw_haemo': raw_haemo,
                    'raw_intensity': raw_intensity,
                    'f_path': path
                }
                
                all_data.append(epoch_data)

    #TODO: raw_haemo throws a weird error when you try to put it into the dataframe
    # dataframe would be better but 
    # all_data_df = pd.DataFrame(all_data)
    # print(all_data_df)
    
    return all_epochs, all_data


def extract_all_amplitudes(all_epochs, tmin=False, tmax=False):
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
                # If tmins and tmax are not specified calculate the event duration from the epoch
                if tmin == False and tmax == False:
                    tmin = subj_data.times[0]
                    tmax = subj_data.times[-1]
                
                value = data.crop(tmin=tmin, tmax=tmax).data * 1.0e6
                # Reshape the data to be a flat numpy array
                value = np.reshape(value, -1)
                
                temporal_measurements.append(value)
    
    temporal_measurements = np.array(temporal_measurements)
    measurement_df = pd.DataFrame(temporal_measurements)

    return measurement_df


def extract_average_amplitudes(all_epochs, tmin=None, tmax=None):
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

    for idx, epoch in enumerate(all_epochs):
        subj_id = 0
        for subj_data in all_epochs[epoch]:
            subj_id += 1
            # can be either "hbo", "hbr", or both
            for chroma in ["hbo", "hbr"]:
                data = deepcopy(subj_data.average(picks=chroma))

                # If tmins and tmax are not specified calculate the event duration from the epoch
                if tmin == None and tmax == None:
                    tmin = subj_data.times[0]
                    tmax = subj_data.times[-1]

                value = data.crop(tmin=tmin, tmax=tmax).data.mean() * 1.0e6

                # Append metadata and extracted feature to dataframe
                this_df = pd.DataFrame(
                    {'ID': subj_id, 'Chroma': chroma, 'Condition': epoch, 'Value': value}, index=[0])
                df = pd.concat([df, this_df], ignore_index=True)

    df['Value'] = pd.to_numeric(df['Value'])  # some Pandas have this as object
    return df


def extract_channel_values(all_epochs, tmin = None, tmax = None):
    """
    Takes in the all_epochs dict and returns a dataframe with the average hemoglobin concentration
    per channel in each condition.
    """
    #TODO: This is a really ugly function. Four nested for loops is pretty bad.
    # We'll have to come back to this to see if we can make it less of fuster cluck
    row_names = {}
    df = pd.DataFrame()
    df_row_number = 0
    for id, event_type in enumerate(all_epochs):
        # Goes through all the epochs generated for the type of task, like "Control"
        for epoch_index, epoch in enumerate(all_epochs[event_type]):
        
            # If no tmin and tmax is given we use the event triggers to delineate the event duration
            if tmin == None and tmax == None:
                tmin = epoch.times[0]
                tmax = epoch.times[-1]

            # Establish were the event occured and crop the rest
            data = epoch.get_data(tmin=tmin, tmax=tmax)
            # test = epoch.to_data_frame()
            if len(epoch.events) >= 1:
                channel_averages = []
                # Loop through each channel and extract the data from that channel
                for channel_index, channel_name in enumerate(epoch.ch_names):
                    # Get the average reading for the channel
                    average = data[:, channel_index, :].mean() * 1.0e6
                    channel_averages.append(average)
                # Create a dataframe to store the channel averages.
                # channel_data[f'{event_type}-{id}'] = channel_averages
                this_df = pd.DataFrame(np.array(channel_averages).reshape(1,-1), columns=epoch.ch_names)        
                # print(this_df)
                df = pd.concat([df, this_df], ignore_index=True)
                # Create a dictionary of row names in order to rename the dataframe indexs
                row_names[df_row_number] = f'{event_type}-{epoch_index + 1}'
                df_row_number += 1
    # Convert indexes to names for ease of use later.
    df = df.rename(index=row_names)
    return df


def find_significant_channels(df_cha):
    mask = df_cha['Significant'] == True
    sig_df = df_cha[mask]
    return sig_df