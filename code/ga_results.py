import numpy as np
import pandas as pd

import glob
import json
from collections import defaultdict

import scipy

import matplotlib.pyplot as plt
import seaborn as sns



def load_ga_results(folder):
    all_files = sorted(glob.glob(folder + "/*.json"))
    
    selected_ppis = []
    selected_corr = []
    score_per_run = defaultdict(list)

    max_corr = 0.0
    max_array = []
    for filename in all_files:
        print(filename)
        with open(filename) as json_file:
            data = json.load(json_file)
            runs = sorted([int(r.lstrip('run')) for r in data.keys()])
            selected_ppis.append(data['run{}'.format(runs[-1])]['ppi_data'])
            selected_corr.append(data['run{}'.format(runs[-1])]['obj_fn'])
            for r in runs:
                corr = data['run{}'.format(r)]['obj_fn']
                score_per_run[r].append(corr)
                if corr > max_corr:
                    max_corr = corr
                    max_array = data['run{}'.format(r)]['ppi_data']
    ga_results = (score_per_run, selected_ppis, selected_corr, max_corr, max_array)                
    return ga_results


def organize_ga_results(ga_results, initial_corr=None):
    score_dict = ga_results[0]
    
    mean_values = []
    std_values = []
    
    if initial_corr is not None:
        mean_values.append(initial_corr)
        std_values.append(0.0)


    for k in sorted(list(score_dict.keys())):
        v = score_dict[k]
        mean_values.append(np.nanmean(v))
        std_values.append(np.nanstd(v))
        
    if initial_corr is not None:
        run_values = list(range(len(mean_values)))
    else:
        run_values = list(range(1, len(mean_values)+1))

    data = np.array(list(zip(run_values, mean_values, std_values)))
    df = pd.DataFrame(data, columns=['X', 'MEAN', 'STD'])
    return df


def plot_run_results(df_list, label_list, filename=None):
    with sns.axes_style("darkgrid"):
        count = 0
        for df, label in zip(df_list, label_list):
            if count == 0:
                ax = plt.plot('X', 'MEAN', data=df, linestyle='', marker='o', markersize=8, 
                             #color='b', 
                             label=f"Mean ({label})"
                            )
            else:
                ax = ax.plot('X', 'MEAN', data=df, linestyle='', marker='o', markersize=8, 
                             #color='b', 
                             label=f"Mean ({label})"
                            )
            mean_precision = df['MEAN'].values
            std_precision = df['STD'].values
            precision_upper = np.minimum(mean_precision + 1.96*std_precision, 1)
            precision_lower = np.maximum(mean_precision - 1.96*std_precision, 0)

            plt.fill_between(df['X'].values, precision_lower, precision_upper, 
                             #color='steelblue', 
                             alpha=.2, 
                             label=f"95% CI ({label})")

            plt.xticks(df['X'].values.tolist())
            plt.tick_params(axis='both', which='major', labelsize=12)
            plt.ylim(0,1)

        plt.xlabel('Iterations of Genetic Algorithm', fontsize=16)
        plt.ylabel('Abs(Spearman score)', fontsize=16)

        plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),
              ncol=1, fancybox=True, shadow=False, fontsize=12)

        if filename is not None:
            plt.savefig(filename,
                        dpi=300,
                        bbox_inches='tight')
    