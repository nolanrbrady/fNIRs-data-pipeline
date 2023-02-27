# Import common libraries
import numpy as np
import pandas as pd
import importlib

# Import MNE processing
import mne

# Import MNE processing
from mne.preprocessing.nirs import optical_density, beer_lambert_law

# Import MNE-NIRS processing
from mne_nirs.statistics import run_glm
from mne_nirs.experimental_design import make_first_level_design_matrix
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




def group_topological_visualisation(df_cha, columns_for_glm_contrast, raw_haemo):
    chromas = ['hbo', 'hbr']
    conditions = len(columns_for_glm_contrast)
    ratios = [1 for i in range(conditions)]

    fig, axes = plt.subplots(nrows=2, ncols=conditions, figsize=(10, 10),
                         gridspec_kw=dict(width_ratios=ratios))
    
    for chroma_id, chroma in enumerate(chromas):
        cmap_color = mpl.cm.Oranges if chroma == 'hbo' else mpl.cm.Blues_r
        vlim = (0, 20) if chroma == 'hbo' else (-10, 0)

        # Cut down the dataframe just to the conditions we are interested in
        ch_summary = df_cha.query(f"Condition in {columns_for_glm_contrast}")
        ch_summary = ch_summary.query(f"Chroma in ['{str(chroma)}']")

        # Run group level model and convert to dataframe
        ch_model = smf.mixedlm("theta ~ -1 + ch_name:Condition:Chroma",
                            ch_summary, groups=ch_summary["ID"]).fit(method='nm')

        ch_model_df = statsmodels_to_results(ch_model)

        for condition_id, condition in enumerate(columns_for_glm_contrast):
            color_bar = True if condition_id + 1 == len(columns_for_glm_contrast) else False
            # Plot the two conditions
            plot_glm_group_topo(raw_haemo.copy().pick(picks=chroma),
                                ch_model_df.query(f"Condition in ['{condition}']"),
                                colorbar=color_bar, axes=axes[chroma_id, condition_id],
                                vlim=vlim, cmap=cmap_color)



def group_plot_visualisation(df_cha, columns_for_glm_contrast, raw_haemo):
    grp_results = df_cha.query("Condition in ['Control', 'Inflam', 'Neutral']")
    grp_results = grp_results.query("Chroma in ['hbo']")
    # print(grp_results)
    cha_model = smf.mixedlm("theta ~ -1 + Condition", grp_results, groups=grp_results["ID"]).fit(method='nm')
    
    df = statsmodels_to_results(cha_model)
    print(df)

    sns.catplot(x="Condition", y="Coef.", hue="Significant", data=df, ci=None, palette="muted", height=4, s=10)


def group_cortical_surface_projection(sub_dir, raw_haemo):
    brain = mne.viz.Brain('fsaverage', subjects_dir=sub_dir, background='w', cortex='0.5')
    brain.add_sensors(raw_haemo.info, trans='fsaverage', fnirs=['channels', 'pairs', 'sources', 'detectors'])
    brain.show_view(azimuth=180, elevation=80, distance=450)