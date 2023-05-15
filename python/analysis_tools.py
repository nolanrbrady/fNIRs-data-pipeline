
# Import MNE-NIRS processing
import mne
from mne_nirs.channels import get_long_channels, get_short_channels
from mne_nirs.datasets import fnirs_motor_group
from mne_nirs.signal_enhancement import enhance_negative_correlation
from mne_nirs.io.snirf import read_snirf_aux_data

# Import MNE processing
from mne import Epochs, events_from_annotations, set_log_level
from mne_nirs.statistics import statsmodels_to_results

# Other Tooling
import pandas as pd
import numpy as np
import os
from collections import defaultdict
from copy import deepcopy
from itertools import compress
import matplotlib.pyplot as plt
import importlib
import statsmodels.formula.api as smf
import statsmodels as sm
from scipy import stats
from statsmodels.stats.weightstats import ztest as ztest


# Local functions
import quality_eval
import dynamic_interval_tools
import custom_timestamp

importlib.reload(quality_eval)
importlib.reload(dynamic_interval_tools)
importlib.reload(custom_timestamp)


def aggregate_epochs(custom_trigger_df, paths, trigger_id, variable_epoch_time, tmin, tmax):
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
    all_evokeds = defaultdict(list)
    # Temporary storage for the items we're using in all_data_df
    all_data = []
    for f_path in paths:

        # Find subject ID from the f_path
        ls = f_path.split('/')
        res = list(filter(lambda a: 'sub' in a, ls))
        sub_id = int(res[0].split('-')[-1])

        epochs, raw_haemo, raw_intensity, path, events, event_dict, aux_df = individual_analysis(
            custom_trigger_df, f_path, trigger_id, variable_epoch_time, tmax, sub_id)

        
        for epoch in epochs:
            for cidx, condition in enumerate(epoch.event_id):
                all_epochs[condition].append(epoch[condition])
                all_evokeds[condition].append(epoch[condition].average())
                epoch_data = {
                    'epoch': epoch,
                    'condition': condition,
                    'raw_haemo': raw_haemo,
                    'raw_intensity': raw_intensity,
                    'f_path': path,
                    'ID': sub_id,
                    'aux_df': aux_df
                }

                all_data.append(epoch_data)
    # TODO: raw_haemo throws a weird error when you try to put it into the dataframe
    # dataframe would be better but
    # all_data_df = pd.DataFrame(all_data)
    # print(all_data_df)

    return all_epochs, all_data, all_evokeds



def individual_analysis(custom_trigger_df, bids_path, trigger_id, variable_epoch_time, tmax, sub_id, custom_trigger=False):
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

    adjusted_path = '/'.join(bids_path.split('/')[:-1])
    content = os.listdir(adjusted_path)

    if len(content) == 2:
       #TODO: concatinate the scans, normalize them and save a consolidated snirf file. 
       print("test", bids_path)

    raw_intensity = mne.io.read_raw_snirf(
        bids_path, verbose=True, preload=False)

    raw_intensity = get_long_channels(raw_intensity, min_dist=0.01)

    if trigger_id:
        
        # Rename the numeric triggers for ease of processing later
        raw_intensity.annotations.rename(trigger_id)

    raw_haemo = quality_eval.signal_preprocessing(raw_intensity)

    # Apply further data cleaning techniques and extract epochs
    raw_haemo = enhance_negative_correlation(raw_haemo)

    # Get Accelerometer data
    aux_df = read_snirf_aux_data(bids_path, raw_haemo)
    

    if custom_trigger_df:
        # TODO: Need to add something here to be able to work in custom triggers
        print('We need to add code to handle custom triggers')
        custom_timestamp.create_custom_events(custom_trigger_df, sub_id)
    else:
        events, event_dict = events_from_annotations(raw_haemo, verbose=False)
        # print(events, event_dict)

    # Logic splits here since there are fundamental differences in how we handle Epoch
    # generation in dynamic intervals instead of block intervals.
    if variable_epoch_time:
        epochs = dynamic_interval_tools.epoch_generation(
            raw_haemo, event_dict, events)
    else:
        # Remove all STOP triggers to hardcode duration to 30 secs per MNE specs
        # TODO: Need a prompt for end triggers present. If NOT present skip this part.
        # events = events[::2]

        epochs = Epochs(raw_haemo, events, event_id=event_dict, tmin=-1, tmax=tmax,
                        reject=dict(hbo=200e-6), reject_by_annotation=True,
                        proj=True, baseline=(None, 0), detrend=0,
                        preload=True, verbose=False)
        # Doing this to ensure that if variable_epoch_time is True or False the format
        # of the data remains the same.
        epochs = [epochs]

    return epochs, raw_haemo, raw_intensity, bids_path, events, event_dict, aux_df


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


def extract_channel_values(all_epochs, tmin=None, tmax=None):
    """
    Takes in the all_epochs dict and returns a dataframe with the average hemoglobin concentration
    per channel in each condition.
    """
    # TODO: This is a really ugly function. Four nested for loops is pretty bad.
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
                this_df = pd.DataFrame(np.array(channel_averages).reshape(
                    1, -1), columns=epoch.ch_names)
                # print(this_df)
                df = pd.concat([df, this_df], ignore_index=True)
                # Create a dictionary of row names in order to rename the dataframe indexs
                row_names[df_row_number] = f'{event_type}-{epoch_index + 1}'
                df_row_number += 1
    # Convert indexes to names for ease of use later.
    df = df.rename(index=row_names)
    return df


def find_significant_channels(ch_model_df, adjust_coef=False):
    keys = list(ch_model_df.keys())
    all_sig_channels = pd.DataFrame()

    for key in keys:
        data = ch_model_df[key]
        mask = data['Significant'] == True
        sig_df = data[mask]
        if adjust_coef:
            sig_df['Coef.'] = sig_df['Coef.']*1e7
        all_sig_channels = all_sig_channels.append(sig_df)

    return all_sig_channels


def create_results_dataframe(df_cha, columns_for_glm_contrast, raw_haemo):
    """
    Returns a dataframe summarizing the significance of all the conditions and channels.
    """
    chromas = ['hbo', 'hbr']

    results = {}

    for chroma in chromas:
        ch_summary = df_cha.query(f"Condition in {columns_for_glm_contrast}")
        ch_summary = ch_summary.query(f"Chroma in ['{chroma}']")

        # Run group level model and convert to dataframe
        ch_model = smf.mixedlm("theta ~ -1 + ch_name:Chroma:Condition",
                               ch_summary, groups=ch_summary["ID"]).fit(method='nm')

        # Here we can use the order argument to ensure the channel name order
        ch_model_df = statsmodels_to_results(ch_model,
                                             order=raw_haemo.copy().pick(
                                                 picks=chroma).ch_names)
        # And make the table prettier
        ch_model_df.reset_index(drop=True, inplace=True)
        ch_model_df = ch_model_df.set_index(['ch_name', 'Condition'])
        results[chroma] = ch_model_df
    return results


def cross_condition_comparison(df, ignore, columns_for_fdr):
    channels = pd.Series(df['ch_name'].unique())
    conditions = pd.Series(df['Condition'].unique())

    compared_channels = pd.DataFrame()

    for channel in channels:
        ch_df = df.loc[df['ch_name'] == channel]

        # Remove the drift columns
        for item in ignore:
            ch_df = ch_df.loc[~ch_df['Condition'].str.contains(item)]

        condition_channels = pd.DataFrame()
        for condition in conditions:
            if (condition not in ignore) and (condition.find('drift') != False) and (condition in columns_for_fdr):
                # Isolate the one condition for the one channel
                condition_df = ch_df.loc[ch_df['Condition'] == condition]

                # Select only the numeric columns
                numeric_columns = condition_df.select_dtypes(
                    include=[float, int])
                non_numeric_columns = condition_df.select_dtypes(exclude=[
                                                                 float, int])

                # Calculate the mean of each numeric column
                column_means = numeric_columns.mean()
                
                # Create a new DataFrame to store the column means in one row
                df_means = pd.DataFrame(column_means).T
                df_numeric = non_numeric_columns.head(1)

                # Isolate the fields we care about
                ch_name = df_numeric['ch_name'].iloc[0]
                chroma = df_numeric['Chroma'].iloc[0]
                condition = df_numeric['Condition'].iloc[0]

                # Add variables to the means
                df_means['ch_name'] = ch_name
                df_means['Chroma'] = chroma
                df_means['Condition'] = condition

                condition_channels = pd.concat([df_means, condition_channels], ignore_index=True)

            # Perform two-sample t-test on the conditions to see if they're different
            if len(condition_channels) > 1:
                cond1 = condition_channels['p_value'][0]
                cond2 = condition_channels['p_value'][1]
                p_vals = [cond1, cond2]

                # Run the T-test to compare the conditions
                print(p_vals)
                result = stats.combine_pvalues(p_vals, method='fisher')
                print(result)
                result_pval = result[1]
                # print("P-val:", result_pval, "Channel:", condition_channels)
                if result_pval <= 1.05:
                    # print("VALID", condition_channels)
                    # fdr_valid_channels = fdr_valid_channels.append(condition_channels)
                    compared_channels = pd.concat([compared_channels, condition_channels], ignore_index=True)

            compared_channels = compared_channels.drop_duplicates()
                    
    return compared_channels
