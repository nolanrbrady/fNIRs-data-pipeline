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

def create_custom_events(csv, sub_id):
    df = pd.read_csv(csv)
    print("Subject ID: ", sub_id)
    # Clean the empty columns
    df = df.filter(regex='^(?!Unnamed).*$', axis=1)
    # Keep only the timestamps from the subject getting processed
    subject_timestamps = df.loc[df['subject_id'] == sub_id]
    print(subject_timestamps)
