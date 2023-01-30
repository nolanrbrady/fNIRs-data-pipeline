import os
import numpy as np
import pandas as pd


def import_data_folder (path, ignored_dirs):
    """
    This function takes in one folder and determines the studies groups and subjects through the fil structure
    """
    print("importing data folder firing: ", path)
    groups = os.listdir(path)
    columns=['group', 'sub_name', 'snirf_path']
    data = []

    # Go through all the groups
    for group in groups:
        if group not in ignored_dirs:
            subs = os.listdir(f'{path}/{group}')

            # Loop through the subjects and extract the snirf files
            for sub in subs:
                if sub not in ignored_dirs:
                    sub_path = f'{path}/{group}/{sub}/nirs'
                    sub_nirs_directory = os.listdir(sub_path)
                    snirf_path = make_snirf_path(sub_nirs_directory, sub_path)
                    data.append([group, sub, snirf_path])

    df = pd.DataFrame(data=data, columns=columns)
    return groups, df




def make_snirf_path(dir,sub_path):
    """
    Simply finds and returns the path of the snirf file in the subjects folder
    """
    for file in dir:
        target = '.snirf'
        if target in file:
            snirf_path = f'{sub_path}/{file}'
            return snirf_path
            