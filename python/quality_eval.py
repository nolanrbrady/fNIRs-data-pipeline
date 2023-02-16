import os
import mne
import numpy as np
from itertools import compress
import matplotlib.pyplot as plt

from mne.preprocessing.nirs import optical_density, temporal_derivative_distribution_repair, scalp_coupling_index
from mne_nirs.visualisation import plot_timechannel_quality_metric
from mne.preprocessing.nirs import beer_lambert_law

def generate_raw_intensity(path):
    raw_intensity = mne.io.read_raw_snirf(path, verbose=True)
    raw_intensity.load_data().resample(4.0, npad="auto")
    return raw_intensity


def evaluate_raw_signal(path):
    raw = generate_raw_intensity(path)
    print(raw.times[-1])
    # Plot the Raw Optical Density
    # raw_od.plot(n_channels=55, duration=4000, show_scrollbars=True, clipping=None)


def evaluate_sci(path):
    raw = generate_raw_intensity(path)
    scan_length = raw.times[-1]
    raw_od = optical_density(raw)
    sci = mne.preprocessing.nirs.scalp_coupling_index(raw_od)

    # Highlights the bad sources and maps a diagram of the montage
    raw_od.info['bads'] = list(compress(raw_od.ch_names, sci < 0.7))
    print(raw_od.info['bads'])

    raw_od.plot_sensors()

    raw_od.copy().plot(n_channels=55, duration=scan_length, show_scrollbars=False,
         clipping=None, scalings={'fnirs_od': 0.2})

    # Plots a histogram of the SCI values
    # fig, ax = plt.subplots()
    # ax.hist(sci)
    # ax.set(xlabel=f'Scalp Coupling Index', ylabel='Count', xlim=[0, 1])
    # plt.show()

def signal_preprocessing(raw_intensity):
    """
    Takes in raw intensity and applies scalp coupling indexing to identify bad channels.
    Once bad channels are removed downsampling occurs, the baseline shift is corrected
    and spike artifacts are removed.

    Returns cleaned raw optical density data.
    """
    # Convert raw signal to optical density
    raw_od = optical_density(raw_intensity)

    # Evaluating the channels for their quality and removing them
    sci = scalp_coupling_index(raw_od, h_freq=1.35, h_trans_bandwidth=0.1)
    raw_od.info["bads"] = list(compress(raw_od.ch_names, sci < 0.5))

    # TODO: Interpolate bads has caused issues in the past but is needed to clean the bad channels.
    # raw_od.interpolate_bads()

    # Downsample and apply signal cleaning techniques
    raw_od.resample(0.8)

    # This approach removes baseline shift and spike artifacts
    raw_od = temporal_derivative_distribution_repair(raw_od)

    # Convert to haemoglobin and filter
    raw_haemo = beer_lambert_law(raw_od, ppf=0.1)
    raw_haemo = raw_haemo.filter(0.02, 0.3,
                                 h_trans_bandwidth=0.1, l_trans_bandwidth=0.01,
                                 verbose=False)
    return raw_haemo

