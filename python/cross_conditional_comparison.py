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
import matplotlib.colors as mcolors
import seaborn as sns

def two_sample_permutation_test(group_data, raw_haemo, columns_for_contrast, contrasts_dict):
    ch_names = raw_haemo.ch_names
    df_con_1_2 = contrasts_dict['contrast_1']
    df_con_2_1 = contrasts_dict['contrast_2']
    contrasts = [df_con_1_2, df_con_2_1]
    conditions = [[0,-1], [-1,0]]
    chromas = ['hbo', 'hbr']
    for idx, contrast in enumerate(contrasts):
        for chroma in chromas:
            cmap = mpl.cm.Oranges if chroma == 'hbo' else mpl.cm.Blues
            # Figure out what the comparison is and generate title
            condition_order = conditions[idx]
            analysis = f'{columns_for_contrast[condition_order[0]]}-{columns_for_contrast[condition_order[-1]]}'
            
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
            con_summary = contrast.query(f"Chroma in ['{str(chroma)}']")

            # Run group level model and convert to dataframe
            con_model = smf.mixedlm("effect ~ -1 + ch_name:Chroma",
                                    con_summary, groups=con_summary["ID"]).fit(method='nm')
        
            contrast_title = f'{analysis} Conditional Difference ({chroma})'
            fig.suptitle(contrast_title)


            con_model_df = statsmodels_to_results(con_model,
                                                order=raw_haemo.copy().pick(
                                                    picks=chroma).ch_names)
            raw = raw_haemo.copy().pick_types(fnirs=True, exclude='bads')
            plot_glm_group_topo(raw.copy().pick(picks=chroma),
                                con_model_df, cmap=cmap, names=ch_names, colorbar=True, axes=axes)

            # ------------------------------------
            # Apply FDR Correction
            # ------------------------------------
            alpha = 0.05
            mask = con_model_df['Significant'] == True
            contrast_df = con_model_df[mask]
            # p_vals = contrast_df['p_value']
            p_vals = contrast_df['P>|z|']
            reject_fdr, pval_fdr = fdr_correction(p_vals, alpha=alpha, method='indep')
            contrast_df['fdr_status'] = reject_fdr

            # Find all the channels that survive FDR
            contrast_df_fdr = contrast_df.loc[(contrast_df['fdr_status'] == True)]

            # NOTE: To see all significant channels across all participants simply use `.drop_duplicates()`
            contrast_df_fdr = contrast_df_fdr.drop_duplicates(subset=['ch_name'], keep='last')

            contrast_df_fdr.to_csv(f'{analysis}_{chroma}_contrast_model_data.csv')
        

