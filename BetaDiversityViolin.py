import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import percentileofscore
import numpy as np


def get_beta_diversity_violin_plot(pickle_file, code, ax=None, excel_writer=None, solo=False):
    combination_names_reference, names, values, actual_values_names, actual_values = pickle_file

    if code == 0:
        names, actual_values_names, values, actual_values = name_converter(names, actual_values_names, values,
                                                                           actual_values)
        if not solo:
            fig, axs = plt.subplots(nrows=4, ncols=7, figsize=(70, 40))
            for combination in range(len(combination_names_reference)):
                axs[int(combination // 7), combination % 7].set_ylabel('Beta Diversity', fontsize=24)
                axs[int(combination // 7), combination % 7].set_ylim(-.1, 1.3)
                axs[int(combination // 7), combination % 7].set_xlim(-0.25, 5.25)
                axs[int(combination // 7), combination % 7].set_xticks([0, 1, 2, 3, 4, 5])
                axs[int(combination // 7), combination % 7].set_yticks([0, .2, .4, .6, .8, 1])
                axs[int(combination // 7), combination % 7].tick_params(axis='y', labelsize=24)
                axs[int(combination // 7), combination % 7].set_xticklabels([], fontsize=24, rotation=45)
                axs[int(combination // 7), combination % 7].set_title(
                    combination_names_reference[combination], fontsize=24)
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
            fig.suptitle('B', x='.05', fontsize=48, fontweight=500)

            ax.set_ylabel('Beta Diversity', fontsize=24)
            ax.set_ylim(-.1, 1.3)
            ax.set_xlim(-0.25, 14.25)
            ax.set_xticks([0, 1, 2, 3, 4, 5])
            ax.set_yticks([0, .2, .4, .6, .8, 1])
            ax.tick_params(axis='y', labelsize=24)
            ax.set_xticklabels([], fontsize=24, rotation=45)
    else:
        if ax is None:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
            fig.suptitle('C', x='.05', fontsize=48, fontweight=500)
        else:
            pass
        ax.set_ylabel('Beta Diversity', fontsize=24)
        ax.set_ylim(-.1, 1.3)
        ax.set_xlim(-0.25, 14.25)
        ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
        ax.set_yticks([0, .2, .4, .6, .8, 1])
        ax.tick_params(axis='y', labelsize=24)
        ax.set_xticklabels([], fontsize=24, rotation=45)

    for combination in range(len(combination_names_reference)):
        plot = True
        if solo and combination != 3:
            plot = False
        if code == 0 and not solo:
            ax = axs[int(combination // 7), combination % 7]
        if plot:
            sns.violinplot(orient='v', x=names[combination], y=values[combination], ax=ax, color='tab:gray')
            sns.swarmplot(orient='v', x=actual_values_names[combination][:], y=actual_values[combination], ax=ax, color='r',
                          size=20)
            percentile_booleans, percentile_booleans_names = get_percentiles(names[combination], values[combination],
                                                                             actual_values_names[combination],
                                                                             actual_values[combination])
            sns.swarmplot(orient='v', x=percentile_booleans_names, y=percentile_booleans, ax=ax, color='k', size=20,
                          marker="*")

    if excel_writer is not None:
        if code == 0:
            if not solo:
                pd.DataFrame(values).transpose().set_axis(combination_names_reference, axis=1).to_excel(excel_writer,
                                                                                                        sheet_name='Fig S3 '
                                                                                                                   'Simulat'
                                                                                                                   'ed Valu'
                                                                                                                   'es')
                pd.DataFrame(actual_values).transpose().set_axis(combination_names_reference, axis=1).to_excel(excel_writer,
                                                                                                        sheet_name='Fig S3 '
                                                                                                                   'Actual '
                                                                                                                   'Values')
        else:
            pd.DataFrame(values).transpose().set_axis(combination_names_reference, axis=1).to_excel(excel_writer,
                                                                                                    sheet_name='Fig 3B '
                                                                                                               'Simulat'
                                                                                                               'ed Valu'
                                                                                                               'es')
            pd.DataFrame(actual_values).transpose().set_axis(combination_names_reference, axis=1).to_excel(excel_writer,
                                                                                                    sheet_name='Fig 3B '
                                                                                                               'Actual '
                                                                                                               'Values')
    if code == 0:
        if solo:
            plt.savefig('/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/'
                        'Fig 5B Beta Diversity Violin Plot.tiff', dpi=300, format='tiff')
        else:
            plt.savefig('/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/'
                        'Fig 3 supplementary Beta Diversity Violin Plots.tiff', dpi=300, format='tiff')
    else:
        pass


def get_percentiles(names, values, actual_names, actual_values, starred=True, percentile_value=5):
    df = pd.DataFrame({'Names': names,
                       'Values': values})
    unique_names = df['Names'].unique()
    grouped_values = [[] for n in unique_names]
    g = df.groupby(['Names'])
    for name in range(len(unique_names)):
        grouped_values[name] += g.get_group(unique_names[name])['Values'].to_list()

    percentile_booleans = []
    percentile_booleans_names = []
    if starred:
        percentiles = [0 for n in range(len(actual_names))]
        for name in range(len(actual_names)):
            index = unique_names.tolist().index(actual_names[name])
            percentiles[name] = percentileofscore(grouped_values[index], actual_values[name])
        for percentile in range(len(percentiles)):
            booleans = [-100, -100]
            if percentiles[percentile] >= 99.0:
                booleans = [1.1, 1.2]
            elif percentiles[percentile] >= 95.0:
                booleans = [1.1, -100]
            percentile_booleans += booleans
            percentile_booleans_names += [actual_names[percentile], actual_names[percentile]]
    else:
        percentile_booleans = []
        percentile_booleans_names = []
        for name in range(len(unique_names)):
            percentile_booleans.append(np.percentile(grouped_values[name], percentile_value))
            percentile_booleans_names.append(unique_names[name])
    return percentile_booleans, percentile_booleans_names


def name_converter(names, actual_values_names, values, actual_values):
    name_converter = {'P79': 'P1',
                      'P84': 'P2',
                      'P73': 'P3',
                      'P55': 'P4',
                      'P52': 'P5',
                      'P76': 'P6'}

    names_new = []
    actual_values_names_new = []

    for comparison in range(len(names)):
        if len(names[comparison]) < 6000:
            continue
        p79_index = names[comparison].index('P79')
        p84_index = names[comparison].index('P84')
        p73_index = names[comparison].index('P73')

        names[comparison] = swapPositions(names[comparison], 0, p79_index)
        names[comparison] = swapPositions(names[comparison], 1, p84_index)
        names[comparison] = swapPositions(names[comparison], 2, p73_index)

        values[comparison] = swapPositions(values[comparison], 0, p79_index)
        values[comparison] = swapPositions(values[comparison], 1, p84_index)
        values[comparison] = swapPositions(values[comparison], 2, p73_index)

        p79_index = actual_values_names[comparison].index('P79')
        p84_index = actual_values_names[comparison].index('P84')
        p73_index = actual_values_names[comparison].index('P73')

        actual_values_names[comparison] = swapPositions(actual_values_names[comparison], 0, p79_index)
        actual_values_names[comparison] = swapPositions(actual_values_names[comparison], 1, p84_index)
        actual_values_names[comparison] = swapPositions(actual_values_names[comparison], 2, p73_index)

        actual_values[comparison] = swapPositions(actual_values[comparison], 0, p79_index)
        actual_values[comparison] = swapPositions(actual_values[comparison], 1, p84_index)
        actual_values[comparison] = swapPositions(actual_values[comparison], 2, p73_index)

    for name_list in names:
        names_new.append([name_converter.get(name) for name in name_list])

    for name_list in actual_values_names:
        actual_values_names_new.append([name_converter.get(name) for name in name_list])

    return names_new, actual_values_names_new, values, actual_values


def swapPositions(l, pos1, pos2):
    l[pos1], l[pos2] = l[pos2], l[pos1]
    return l