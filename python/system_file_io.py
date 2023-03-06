import os
import numpy as np
import pandas as pd


def import_data_folder (path, ignored_dirs):
    """
    This function takes in one folder and determines the studies groups and subjects through the fil structure

    File Structure Expected:
    - Data
        - Group 1
            - Sub-1
                - nirs
                    - something.snirf
            - Sub-2
                - nirs
                    - something.snirf
            - Sub-3
                - nirs
                    - something.snirf
        - Group-2
            - Sub-4
                - nirs
                    - something.snirf
            - Sub-5
                - nirs
                    - something.snirf
            - Sub-6
                - nirs
                    - something.snirf
    """
    groups = os.listdir(path)
    groups = list(filter(lambda group: group not in ignored_dirs, groups))
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