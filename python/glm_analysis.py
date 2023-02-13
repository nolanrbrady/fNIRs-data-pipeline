# Import common libraries
import numpy as np
import pandas as pd

# Import MNE processing
from mne.preprocessing.nirs import optical_density, beer_lambert_law

# Import MNE-NIRS processing
from mne_nirs.statistics import run_glm
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne import Epochs, events_from_annotations
from mne_nirs.statistics import statsmodels_to_results
from mne_nirs.channels import get_short_channels, get_long_channels
from mne_nirs.channels import picks_pair_to_idx
from mne_nirs.visualisation import plot_glm_group_topo
from mne_nirs.datasets import fnirs_motor_group
from mne_nirs.visualisation import plot_glm_surface_projection
from mne_nirs.io.fold import fold_channel_specificity

# Import MNE-BIDS processing
from mne_bids import BIDSPath, read_raw_bids, get_entity_vals

# Import StatsModels
import statsmodels.formula.api as smf

# Import Plotting Library
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns


def create_design_matrix(all_data):
    updated_data = []
    for data in all_data:
        epoch, condition, raw_haemo = data.values()
        print(raw_haemo)
        events, event_dict = events_from_annotations(raw_haemo, verbose=False)
        
        for event in events:
            # Dynamically establish the task length
            prev_event_time = events[-2][0]
            current_event_time = events[-1][0]
            task_len = current_event_time - prev_event_time
        print(task_len)
        design_matrix = make_first_level_design_matrix(raw_haemo, stim_dur=task_len)
        data['design_matrix'] = design_matrix
        updated_data.append(data)
        print(data)
    
    return updated_data


def create_glm_df(glm_data, columns_for_contrast=None):
    df_cha = pd.DataFrame() # Stores channel level results
    df_con = pd.DataFrame() # Stores channel level contrast results
    for data in glm_data:
        raw_haemo = data['raw_haemo']
        design_matrix = data['design_matrix']

        # Convert GLM into a Dataframe
        glm_est = run_glm(raw_haemo, design_matrix)
        cha = glm_est.to_dataframe()

        if columns_for_contrast:
            # Define the GLM contrast that is to be evaluated
            contrast_matrix = np.eye(design_matrix.shape[1])
            
            basic_conts = dict([(column, contrast_matrix[i])for i, column in enumerate(design_matrix.columns)])
            column_1 = columns_for_contrast[0]
            column_2 = columns_for_contrast[-1]
            contrast_LvR = basic_conts[column_1] - basic_conts[column_2]

            # Compute defined contrast
            contrast = glm_est.compute_contrast(contrast_LvR)
            con = contrast.to_dataframe()

            df_con = pd.concat([df_con, con], ignore_index=True)

        

        df_cha = pd.concat([df_cha, cha], ignore_index=True)
        
    return df_cha, df_con