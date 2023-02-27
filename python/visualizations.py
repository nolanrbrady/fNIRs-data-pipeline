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
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10),
                         gridspec_kw=dict(width_ratios=[1, 1]))

    # Cut down the dataframe just to the conditions we are interested in
    ch_summary = df_cha.query("Condition in ['Control', 'Neutral', 'Inflam']")
    ch_summary = ch_summary.query("Chroma in ['hbo']")

    # Run group level model and convert to dataframe
    ch_model = smf.mixedlm("theta ~ -1 + ch_name:Condition:Chroma",
                        ch_summary, groups=ch_summary["ID"]).fit(method='nm')

    ch_model_df = statsmodels_to_results(ch_model)

    # Plot the two conditions
    plot_glm_group_topo(raw_haemo.copy().pick(picks="hbo"),
                        ch_model_df.query("Condition in ['Control']"),
                        colorbar=False, axes=axes[0, 0],
                        vlim=(0, 20), cmap=mpl.cm.Oranges)

    plot_glm_group_topo(raw_haemo.copy().pick(picks="hbo"),
                        ch_model_df.query("Condition in ['Inflam']"),
                        colorbar=True, axes=axes[0, 1],
                        vlim=(0, 20), cmap=mpl.cm.Oranges)

    # Cut down the dataframe just to the conditions we are interested in
    ch_summary = df_cha.query("Condition in ['Control', 'Inflam']")
    ch_summary = ch_summary.query("Chroma in ['hbr']")

    # Run group level model and convert to dataframe
    ch_model = smf.mixedlm("theta ~ -1 + ch_name:Condition:Chroma",
                        ch_summary, groups=ch_summary["ID"]).fit(method='nm')
    ch_model_df = statsmodels_to_results(ch_model)

    # Plot the two conditions
    plot_glm_group_topo(raw_haemo.copy().pick(picks="hbr"),
                        ch_model_df.query("Condition in ['Control']"),
                        colorbar=False, axes=axes[1, 0],
                        vlim=(-10, 0), cmap=mpl.cm.Blues_r)
    plot_glm_group_topo(raw_haemo.copy().pick(picks="hbr"),
                        ch_model_df.query("Condition in ['Inflam']"),
                        colorbar=True, axes=axes[1, 1],
                        vlim=(-10, 0), cmap=mpl.cm.Blues_r)


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