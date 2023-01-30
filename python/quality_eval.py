import os
import mne
import numpy as np
from itertools import compress
import matplotlib.pyplot as plt

from mne.preprocessing.nirs import optical_density
from mne_nirs.preprocessing import peak_power, scalp_coupling_index_windowed
from mne_nirs.visualisation import plot_timechannel_quality_metric

def generate_raw_intensity(path):
    raw_intensity = mne.io.read_raw_snirf(path, verbose=True)
    raw_intensity.load_data().resample(4.0, npad="auto")
    return raw_intensity


def evaluate_raw_signal(path):
    raw = generate_raw_intensity(path)

    raw_od = optical_density(raw)

    # Plot the Raw Optical Density
    raw_od.plot(n_channels=55, duration=4000, show_scrollbars=True, clipping=None)


def evaluate_sci(path):
    raw = generate_raw_intensity(path)
    raw_od = optical_density(raw)
    sci = mne.preprocessing.nirs.scalp_coupling_index(raw_od)

    # Highlights the bad sources and maps a diagram of the montage
    raw_od.info['bads'] = list(compress(raw_od.ch_names, sci < 0.7))
    print(raw_od.info['bads'])

    raw_od.plot_sensors()

    raw_od.copy().plot(n_channels=55, duration=40000, show_scrollbars=False,
         clipping=None, scalings={'fnirs_od': 0.2})

    # Plots a histogram of the SCI values
    # fig, ax = plt.subplots()
    # ax.hist(sci)
    # ax.set(xlabel=f'Scalp Coupling Index', ylabel='Count', xlim=[0, 1])
    # plt.show()

