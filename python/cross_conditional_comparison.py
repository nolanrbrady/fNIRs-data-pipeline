# Import common libraries
import numpy as np
import pandas as pd

# Import MNE processing
from mne.preprocessing.nirs import optical_density, beer_lambert_law

# Import MNE-NIRS processing
from mne_nirs.statistics import run_glm
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.statistics import statsmodels_to_results
from mne.stats import fdr_correction
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

def two_sample_permutation_test(group_data, raw_haemo, columns_for_contrast, contrasts_dict):
    ch_names = raw_haemo.ch_names
    column_1 = columns_for_contrast[0]
    column_2 = columns_for_contrast[-1]
    df_con_1_2 = contrasts_dict['contrast_1']
    df_con_2_1 = contrasts_dict['contrast_2']
    contrasts = [df_con_1_2, df_con_2_1]
    conditions = [column_1, column_2]
    # ------------------------------------
    # Contrast Condition 1 - Condition 2
    # ------------------------------------
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    con_1_2_summary = df_con_1_2.query("Chroma in ['hbo']")

    # Run group level model and convert to dataframe
    con_model = smf.mixedlm("effect ~ -1 + ch_name:Chroma",
                            con_1_2_summary, groups=con_1_2_summary["ID"]).fit(method='nm')
    
    analysis = f'{column_1}-{column_2}'
    # NOTE: This is based off of the code in the `create_glm_df()` function
    contrast_title = f'{analysis} Conditional Difference'
    fig.suptitle(contrast_title)


    con_model_df = statsmodels_to_results(con_model,
                                        order=raw_haemo.copy().pick(
                                            picks="hbo").ch_names)

    plot_glm_group_topo(raw_haemo.copy().pick(picks="hbo"),
                        con_model_df, names=ch_names, colorbar=True, axes=axes)
    
    # ------------------------------------
    # Contrast Condition 2 - Condition 1
    # ------------------------------------
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    con_2_1_summary = df_con_2_1.query("Chroma in ['hbo']")

    # Run group level model and convert to dataframe
    con_model = smf.mixedlm("effect ~ -1 + ch_name:Chroma",
                            con_2_1_summary, groups=con_2_1_summary["ID"]).fit(method='nm')
    
    analysis = f'{column_2}-{column_1}'
    # NOTE: This is based off of the code in the `create_glm_df()` function
    contrast_title = f'{analysis} Conditional Difference'
    fig.suptitle(contrast_title)


    con_model_df = statsmodels_to_results(con_model,
                                        order=raw_haemo.copy().pick(
                                            picks="hbo").ch_names)

    plot_glm_group_topo(raw_haemo.copy().pick(picks="hbo"),
                        con_model_df, names=ch_names, colorbar=True, axes=axes)


    # ------------------------------------
    # Apply FDR Correction
    # ------------------------------------
    alpha = 0.05

    for idx, constrast in enumerate(contrasts):
        
        mask = constrast['Significant'] == True
        contrast_df = constrast[mask]
        p_vals = contrast_df['p_value']
        reject_fdr, pval_fdr = fdr_correction(p_vals, alpha=alpha, method='indep')
        contrast_df['fdr_status'] = reject_fdr

        # Find all the channels that survive FDR
        contrast_df_fdr = contrast_df.loc[(contrast_df['fdr_status'] == True)]

        # NOTE: To see all significant channels across all participants simply use `.drop_duplicates()`
        contrast_df_fdr = contrast_df_fdr.drop_duplicates(subset=['ch_name'], keep='last')

        # print(contrast_df_fdr)
        contrast_df_fdr.to_csv(f'{analysis}_contrast_model_data.csv')
    

