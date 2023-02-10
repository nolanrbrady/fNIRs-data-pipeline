
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
from copy import deepcopy

# Other Tooling
import pandas as pd
import numpy as np



def epoch_generation(raw_haemo, event_dict, events):
    all_epochs = []
    adjusted_dict = { 1: {'Practice': 1},
        2: {'Neutral': 2},
        3: {'Inflammatory': 3},
        4: {'Control': 4}}

    for index, event in enumerate(events):
        total_len = len(events) - 1

        if index == total_len:
            event_type = events[- 1][2]
            prev_event_time = events[-2][0]
            current_event_time = events[-1][0]
            task_len = current_event_time - prev_event_time
            event_id = adjusted_dict[event_type]
            
            # Create a temporary events array for the epoch creation
            # Only use the start of the event for the event creation
            # prev_event_trigger is time stamp for the start trigger
            current_event = [[prev_event_time, 0, event_type]]
            

            epochs = Epochs(raw_haemo, events = current_event, event_id=event_id, 
                        tmin=-1, tmax=task_len,
                        reject=dict(hbo=200e-6), reject_by_annotation=True,
                        proj=True, baseline=(None, 0), detrend=0,
                        preload=True, verbose=False) 

            all_epochs.append(epochs) 
            # print(adjusted_dict[event_type])
        
        if index % 2 == 0 and index != 0:
            # Find the event_id to be able to match it with the event_dict
            event_type = events[index - 1][2]
            
            # Identify the first and second triggers time stamps
            prev_event_time = events[index - 2][0]
            current_event_time = events[index - 1][0]

            # Finds the length between two "working" period triggers
            task_len = current_event_time - prev_event_time

            event_id = adjusted_dict[event_type]

            # Create a temporary events array for the epoch creation
            # Only use the start of the event for the event creation
            # prev_event_trigger is time stamp for the start trigger
            current_event = [[prev_event_time, 0, event_type]]
            
            epochs = Epochs(raw_haemo, events = current_event, event_id=event_id, 
                        tmin=-1, tmax=task_len,
                        reject=dict(hbo=200e-6), reject_by_annotation=True,
                        proj=True, baseline=(None, 0), detrend=0,
                        preload=True, verbose=False)        

            all_epochs.append(epochs)
            # print(adjusted_dict[event_type])

        print(all_epochs)

    return all_epochs
