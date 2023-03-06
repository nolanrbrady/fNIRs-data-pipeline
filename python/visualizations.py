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
from mne.viz import plot_compare_evokeds
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
from pprint import pprint




def group_topological_visualisation(df_cha, columns_for_glm_contrast, raw_haemo, group):
    chromas = ['hbo', 'hbr']
    conditions = len(columns_for_glm_contrast)
    ratios = [1 for i in range(conditions)]

    fig, axes = plt.subplots(nrows=2, ncols=conditions, figsize=(10, 10),
                         gridspec_kw=dict(width_ratios=ratios))
    
    fig.suptitle(group.capitalize(), fontsize=20)
    
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

        ch_model_df['Coef.'] = ch_model_df['Coef.'] * 1e7 

        for condition_id, condition in enumerate(columns_for_glm_contrast):
            color_bar = True if condition_id + 1 == len(columns_for_glm_contrast) else False
            # Plot the condition
            plot_glm_group_topo(raw_haemo.copy().pick(picks=chroma),
                                ch_model_df.query(f"Condition in ['{condition}']"),
                                colorbar=color_bar, axes=axes[chroma_id, condition_id],
                                vlim=vlim, cmap=cmap_color, threshold=True)
    plt.show()



def group_plot_visualisation(df_cha, columns_for_glm_contrast, raw_haemo):
    grp_results = df_cha.query("Condition in ['Control', 'Inflam', 'Neutral']")
    grp_results = grp_results.query("Chroma in ['hbo']")
    # print(grp_results)
    cha_model = smf.mixedlm("theta ~ -1 + Condition", grp_results, groups=grp_results["ID"]).fit(method='nm')
    
    df = statsmodels_to_results(cha_model)
    print(df)

    sns.catplot(x="Condition", y="Coef.", hue="Significant", data=df, ci=None, palette="muted", height=4, s=10)


def group_cortical_surface_projection(df_cha, columns_for_glm_constrast, raw_haemo):
    #TODO: This doesn't work. We need to fix it.
    cha_summary = df_cha.query("Chroma in ['hbo']")

    # Run group level model and convert to dataframe
    con_model = smf.mixedlm("theta ~ -1 + ch_name:Chroma",
                            cha_summary, groups=cha_summary["ID"]).fit(method='nm')
    con_model_df = statsmodels_to_results(con_model,
                                        order=raw_haemo.copy().pick(
                                            picks="hbo").ch_names)
    # Generate brain figure from data
    clim = dict(kind='value', pos_lims=(0, 8, 11))
    brain = plot_glm_surface_projection(raw_haemo.copy().pick("hbo"),
                                        con_model_df, clim=clim, view='dorsal',
                                        colorbar=True, size=(800, 700))
    
    brain.add_text(0.05, 0.95, "Left-Right", 'title', font_size=16, color='k')

    # Run model code as above
    clim = dict(kind='value', pos_lims=(0, 11.5, 17))
    for idx, cond in enumerate(columns_for_glm_constrast):

        # Run same model as explained in the sections above
        ch_summary = df_cha.query("Condition in [@cond]")
        ch_summary = ch_summary.query("Chroma in ['hbo']")
        ch_model = smf.mixedlm("theta ~ -1 + ch_name", ch_summary,
                            groups=ch_summary["ID"]).fit(method='nm')
        model_df = statsmodels_to_results(ch_model, order=raw_haemo.copy().pick("hbo").ch_names)

        # Generate brain figure from data
        brain = plot_glm_surface_projection(raw_haemo.copy().pick("hbo"),
                                            model_df, clim=clim, view='dorsal',
                                            colorbar=True, size=(800, 700))
        brain.add_text(0.05, 0.95, cond, 'title', font_size=16, color='k')


def plot_waveform_analysis(all_evokeds, variable_epoch_time):
    if variable_epoch_time:
        return RuntimeError('Waveform Analysis is not possible with variable length tasks')
    else:
        # Specify the figure size and limits per chromophore
        fig, axes = plt.subplots(nrows=1, ncols=len(all_evokeds), figsize=(17, 5))
        lims = dict(hbo=[-5, 12], hbr=[-5, 12])

        for (pick, color) in zip(['hbo', 'hbr'], ['r', 'b']):
            for idx, evoked in enumerate(all_evokeds):
                print(evoked, len(all_evokeds[evoked]))
                plot_compare_evokeds({evoked: all_evokeds[evoked]}, combine='mean',
                                    picks=pick, axes=axes[idx], show=False,
                                    colors=[color], legend=False, ylim=lims, ci=0.95,
                                    show_sensors=idx == 2)
                axes[idx].set_title('{}'.format(evoked))
        axes[0].legend(["Oxyhaemoglobin", "Deoxyhaemoglobin"])
    plt.show()