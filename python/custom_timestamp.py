import numpy as np
import mne
from mne_nirs.channels import get_long_channels
from mne_nirs.channels import picks_pair_to_idx
from mne_nirs.datasets import fnirs_motor_group
from mne.viz import plot_compare_evokeds
from mne.preprocessing.nirs import optical_density, temporal_derivative_distribution_repair, scalp_coupling_index
from mne_nirs.visualisation import plot_timechannel_quality_metric
from mne.preprocessing.nirs import beer_lambert_law
from mne import Epochs, events_from_annotations, set_log_level
from mne.filter import resample
from mne_nirs.signal_enhancement import enhance_negative_correlation
import matplotlib.pyplot as plt
from itertools import compress
from datetime import datetime
import os
import time
import pandas as pd
import numpy as np

def create_custom_events(csv, sub_id, raw_haemo):
    df = pd.read_csv(csv)
    start_time = str(raw_haemo.info['meas_date']).split(' ')[-1]
    start_time = start_time.split('+')[0]
    scan_start_timestamp = datetime.strptime(start_time, "%H:%M:%S")
    sfreq = raw_haemo.info['sfreq']
    # print("Subject ID: ", sub_id, scan_start_timestamp, sfreq)
    # Clean the empty columns
    df = df.filter(regex='^(?!Unnamed).*$', axis=1)
    # Keep only the timestamps from the subject getting processed
    subject_timestamps = df.loc[df['subject_id'] == sub_id]
    # print(subject_timestamps)

    # Adjust the timestamps to be number of scenes from start instead of timestamps
    def convert_timestamp_to_sample(entry):
        trigger_timestamp = datetime.strptime(entry, "%H:%M:%S")
        delta = trigger_timestamp - scan_start_timestamp

        # convert timestamps to seconds then to samples
        hours, minutes, seconds = str(delta).split(':')
        sample = round(((int(hours) * 60 * 60) + (int(minutes) * 60) + int(seconds)) / sfreq)
        return sample
    
    subject_timestamps['start'] = subject_timestamps['start'].apply(convert_timestamp_to_sample)
    subject_timestamps['end'] = subject_timestamps['end'].apply(convert_timestamp_to_sample)

    print(subject_timestamps)
