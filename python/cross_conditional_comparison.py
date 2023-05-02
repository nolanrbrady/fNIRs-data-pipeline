import numpy as np
from scipy import stats as stats

import mne
from mne import spatial_src_adjacency
from mne.channels import find_ch_adjacency
from mne.stats import spatio_temporal_cluster_test, summarize_clusters_stc
from mne import Epochs, events_from_annotations
from mne.datasets import sample


def two_sample_permutation_test(epochs, raw_data, tmin, tmax):
    # Try Passing in Raw intead of Epoch Data
    for raw in raw_data:
        ch_names = raw.ch_names
        events, event_dict = events_from_annotations(raw, verbose=False)

        epoch = mne.Epochs(raw, events, event_dict, tmin, tmax, picks=None,
                    baseline=None, reject=None, preload=True)
        
        epoch.pick_channels(ch_names)
        epoch.equalize_event_counts(event_dict)

        # for channel in ch_names:
        #     pick = []
        #     pick.append(channel)
        #     print(pick)
        #     epoch.pick_channels(channel)
        #     epoch.equalize_event_counts(event_dict)