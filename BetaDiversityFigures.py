import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import SpeciesData as sd
import IndividualData as id
import seaborn as sns
import HeatmapLabeler as hl
from itertools import combinations
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import LabelEncoder
from scipy.stats import ttest_rel, ttest_ind
from scipy.stats import percentileofscore
from scipy import stats
import pickle
from collections import Counter
import pickle


pig_lung_day3_locations = {'L.Apical': 0,
                           'R.Apical': 1,
                           'L.Cardiac': 2,
                           'R.Cardiac': 3,
                           'L.Diaphr': 4,
                           'R.Diaphr': 5,
                           'Intermed': 6}
pig_lung_day5_locations = {'L.Apical': 0,
                           'R.Apical': 1,
                           'L.Diaphr': 2,
                           'R.Diaphr': 3,
                           'L.Cardiac': 4,
                           'R.Cardiac': 5,
                           'Intermed': 6}
species_locations = {'Pig Lung Day 3': 0,
                     'Pig Lung Day 5': 1}
individual_names_locations = {'P79': 0, 'P84': 1, 'P73': 2, 'P55': 3, 'P52': 4, 'P76': 5}


def get_beta_diversity_figure(ds):
    colors = {'0': 'tab:brown', '1': 'tab:olive', '2': 'tab:gray'}
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20, 20))

    for individual in range(ds.number_of_individuals):
        ax1.scatter(0, ds.individuals[individual].beta_diversity_0, s=96, c=colors.get(f'{individual}'), alpha=0.8,
                    label=ds.individual_names[individual])
        ax1.scatter(1, ds.individuals[individual].beta_diversity_1, s=96, c=colors.get(f'{individual}'), alpha=0.8)
        ax1.scatter(2, ds.individuals[individual].beta_diversity_infinity, s=96, c=colors.get(f'{individual}'),
                    alpha=0.8)

        ax2.scatter(0, ds.individuals[individual].Meff_0_standardized, s=96, c=colors.get(f'{individual}'), alpha=0.8)
        ax2.scatter(1, ds.individuals[individual].Meff_1_standardized, s=96, c=colors.get(f'{individual}'), alpha=0.8)
        ax2.scatter(2, ds.individuals[individual].Meff_infinity_standardized, s=96, c=colors.get(f'{individual}'),
                    alpha=0.8)

    ax1.set_xlim(-0.05, 2.05)
    ax1.set_xticks([0, 1, 2])
    ax1.set_xticklabels(['q = 0', 'q = 1', 'q = ∞'], fontsize=30)
    ax1.set_ylabel('Diversity', fontsize=30)
    ax1.legend(fontsize=30)
    ax1.set_title('A', loc='left', fontsize=48, fontweight=500)
    ax1.tick_params(axis='y', labelsize=30)

    ax2.set_xlim(-0.05, 2.05)
    ax2.set_xticks([0, 1, 2])
    ax2.set_xticklabels(['q = 0', 'q = 1', 'q = ∞'], fontsize=30)
    ax2.set_ylim(-0.05, 1.05)
    ax2.set_ylabel('Diversity', fontsize=30)
    ax2.set_title('B', loc='left', fontsize=48, fontweight=500)
    ax2.tick_params(axis='y', labelsize=30)

    plt.savefig(f'/Users/maxbagga/Desktop/Box/Mallard Phylogeny Data/{ds.species_name} Results/{ds.species_name} '
                f'Beta Diversity Figure.png')
    plt.show()


def get_beta_diversity_matrix(ds, is_jaccard_index=False, is_simulation=False, is_ttest=False, solo=True,
                              species_label=None, ncols=None, nrows=None, suptitle=None, code=0):
    """Creates heatmap matrix for the three beta diversities for each individual in dataset"""
    if not is_jaccard_index:
        Meff_0_values = []
        Meff_1_values = []
        Meff_infinity_values = []
    else:
        jaccard_values = []
    names = []
    if ncols is not None:
        if code != 2:
            fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(20, 5))
            fig.suptitle(suptitle, x='.05', fontsize=48, fontweight=500)
        else:
            fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(20, 5), constrained_layout=True)
    ax_count = 0
    for individual in ds.individuals:
        names.append(individual.name)
        if not is_jaccard_index:
            Meff_0_values_individual = []
            Meff_1_values_individual = []
            Meff_infinity_values_individual = []
        else:
            jaccard_values_individual = []
        if not is_jaccard_index:
            Meff_0_standardized = np.zeros((ds.max_len_data_points, ds.max_len_data_points))
            Meff_1_standardized = np.zeros((ds.max_len_data_points, ds.max_len_data_points))
            Meff_infinity_standardized = np.zeros((ds.max_len_data_points, ds.max_len_data_points))
        else:
            jaccard_index = np.zeros((ds.max_len_data_points, ds.max_len_data_points))
        for combination in list(combinations(range(individual.number_of_days), 2)):
            x = combination[0]
            y = combination[1]
            individual_data_temp = id.Individual(individual.wb, individual.sheet, individual.name,
                                                 [individual.x_dimension_starts[x], individual.x_dimension_starts[y]],
                                                 [individual.x_dimension_finishes[x],
                                                  individual.x_dimension_finishes[y]],
                                                 [individual.y_dimension_starts[x], individual.y_dimension_starts[y]],
                                                 [individual.y_dimension_finishes[x],
                                                  individual.y_dimension_finishes[y]],
                                                 [individual.data_points[x], individual.data_points[y]],
                                                 is_simulation=is_simulation)
            if not is_jaccard_index:
                Meff_0_standardized[x][y] = individual_data_temp.Meff_0_standardized
                Meff_0_standardized[y][x] = individual_data_temp.Meff_0_standardized
                Meff_1_standardized[x][y] = individual_data_temp.Meff_1_standardized
                Meff_1_standardized[y][x] = individual_data_temp.Meff_1_standardized
                Meff_infinity_standardized[x][y] = individual_data_temp.Meff_infinity_standardized
                Meff_infinity_standardized[y][x] = individual_data_temp.Meff_infinity_standardized
                Meff_0_values_individual.append(individual_data_temp.Meff_0_standardized)
                Meff_1_values_individual.append(individual_data_temp.Meff_1_standardized)
                Meff_infinity_values_individual.append(individual_data_temp.Meff_infinity_standardized)
            else:
                jaccard_index[x][y] = individual_data_temp.jaccard_index
                jaccard_index[y][x] = individual_data_temp.jaccard_index
                jaccard_values_individual.append(individual_data_temp.jaccard_index)
        if not is_jaccard_index:
            Meff_0_values.append(Meff_0_values_individual)
            Meff_1_values.append(Meff_1_values_individual)
            Meff_infinity_values.append(Meff_infinity_values_individual)
        else:
            jaccard_values.append(jaccard_values_individual)

        if not is_jaccard_index and not is_ttest:
            if solo:
                df_Meff_0_standardized = pd.DataFrame(Meff_0_standardized, index=individual.data_points,
                                                      columns=individual.data_points)

                mask = np.zeros_like(df_Meff_0_standardized)
                mask[np.triu_indices_from(mask)] = True
                if nrows == 1:
                    ax = axs[ax_count]
                else:
                    ax = axs[ax_count // ncols, ax_count % ncols]
                with sns.axes_style("white"):
                    if ax_count != ds.number_of_individuals - 1:
                        sns.heatmap(df_Meff_0_standardized, vmin=0, vmax=1, ax=ax,
                                    cbar_kws={'label': 'Beta Diversity'}, yticklabels=False, xticklabels=False,
                                    mask=mask, cmap='crest', cbar=False)
                    else:
                        sns.heatmap(df_Meff_0_standardized, vmin=0, vmax=1, ax=ax,
                                    cbar_kws={'label': 'Beta Diversity'}, yticklabels=False, xticklabels=False,
                                    mask=mask, cmap='crest', cbar=True)
                if ax_count == 0:
                    if code == 0:
                        ax.plot([], [], ' ', label="1. L.Apical\n2. R.Apical\n3. L.Cardiac\n4. R.Cardiac\n5. "
                                                   "L.Diaphr\n6. R.Diaphr\n7. Intermed\n8. Nasal")
                        ax.legend(loc="upper right", handlelength=0, handletextpad=0)
                    elif code == 1:
                        ax.plot([], [], ' ', label="1. L.Apical\n2. R.Apical\n3. L.Cardiac\n4. R.Cardiac\n5. "
                                                   "L.Diaphr\n6. R.Diaphr\n7. Intermed")
                        ax.legend(loc="upper right", handlelength=0, handletextpad=0)
                if code != 2:
                    ax.set_title(ds.individuals[ax_count].name, fontsize=30)
                else:
                    ax.set_title(ds.individuals[ax_count].name, fontsize=30, loc='left')
                ax_count += 1
            else:
                fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(30, 24/3))

                df_Meff_0_standardized = pd.DataFrame(Meff_0_standardized, index=individual.data_points,
                                                      columns=individual.data_points)
                df_Meff_1_standardized = pd.DataFrame(Meff_1_standardized, index=individual.data_points,
                                                      columns=individual.data_points)
                df_Meff_infinity_standardized = pd.DataFrame(Meff_infinity_standardized, index=individual.data_points,
                                                             columns=individual.data_points)

                sns.heatmap(df_Meff_0_standardized, vmin=0, vmax=1, ax=ax1, cbar_kws={'label': 'Beta Diversity'},
                            annot=True)
                sns.heatmap(df_Meff_1_standardized, vmin=0, vmax=1, ax=ax2, cbar_kws={'label': 'Beta Diversity'},
                            annot=True)
                sns.heatmap(df_Meff_infinity_standardized, vmin=0, vmax=1, ax=ax3, cbar_kws={'label': 'Beta Diversity'},
                            annot=True)

                ax1.set_title('q = 0', fontsize=30)
                ax2.set_title('q = 1', fontsize=30)
                ax3.set_title('q = ∞', fontsize=30)

                if not is_simulation:
                    plt.savefig(f'/Users/maxbagga/Desktop/Box/Mallard Phylogeny Data/{ds.species_name} Results/Beta Diversity/'
                                f'{ds.species_name} Beta Diversity Matrix {individual.name}.png')
                else:
                    plt.savefig(f'/Users/maxbagga/Desktop/Box/Mallard Phylogeny Data/{ds.species_name} Results/Beta Diversity/'
                                f'{ds.species_name} Beta Diversity Matrix {individual.name} Simulated.png')
        elif is_jaccard_index and not is_ttest:
            fig, ax1 = plt.subplots(ncols=1, figsize=(10, 24 / 3))

            df_jaccard_index = pd.DataFrame(jaccard_index, index=individual.data_points, columns=individual.data_points)

            sns.heatmap(df_jaccard_index, vmin=0, vmax=1, ax=ax1, cbar_kws={'label': 'Jaccard Index'},
                        annot=True)

            ax1.set_title('Jaccard Index', fontsize=30)

            if not is_simulation:
                plt.savefig(f'/Users/maxbagga/Desktop/Box/Mallard Phylogeny Data/{ds.species_name} Results/Jaccard Index/'
                            f'{ds.species_name} Jaccard Index Matrix {individual.name}.png')
            else:
                plt.savefig(f'/Users/maxbagga/Desktop/Box/Mallard Phylogeny Data/{ds.species_name} Results/Jaccard Index/'
                            f'{ds.species_name} Jaccard Index Matrix {individual.name} Simulated.png')

        # plt.show()
    if ncols is not None:
        plt.savefig(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/'
                    f'Fig 2 supplementary {species_label}.tiff', dpi=300, format='tiff')
    if not is_jaccard_index:
        return Meff_0_values, Meff_1_values, Meff_infinity_values, names
    else:
        return jaccard_values


def get_beta_diversity_matrix_across_indiviuals(ds, is_across_ds=False, is_jaccard_index=False, is_pca=False,
                                                is_simulation=False, specific_lobe=None, solo=True, species_label=None,
                                                ax=None, is_cube=False, excel_writer=None):
    """Creates heatmap matrix for the three beta diversities across individuals in dataset"""
    individuals_number_of_days = []
    names = []
    Meff_0_values = []
    Meff_1_values = []
    Meff_infinity_values = []
    jaccard_values = []
    if not is_across_ds:
        for individual in ds.individuals:
            individuals_number_of_days.append(individual.number_of_days)
            names += [f'{individual.name} {day.day}' for day in individual.days]
    else:
        for individual in range(ds.number_of_individuals):
            individuals_number_of_days.append(ds.individuals[individual].number_of_days)
            if individual < ds.ds_split_point:
                names += [f'{ds.species_name[0]} {ds.individuals[individual].name} {day.day}' for day in
                          ds.individuals[individual].days]
            else:
                names += [f'{ds.species_name[1]} {ds.individuals[individual].name} {day.day}' for day in
                          ds.individuals[individual].days]
    population_size = sum(individuals_number_of_days)
    individuals = []
    for individual in range(len(individuals_number_of_days)):
        if individual == 0:
            individuals += [range(individuals_number_of_days[individual])]
        else:
            individuals += [range(individuals[individual - 1][-1] + 1, individuals[individual - 1][-1] +
                                  individuals_number_of_days[individual] + 1)]

    if not is_pca:
        if not is_jaccard_index:
            Meff_0_standardized = np.zeros((population_size, population_size))
            Meff_1_standardized = np.zeros((population_size, population_size))
            Meff_infinity_standardized = np.zeros((population_size, population_size))
        else:
            jaccard_index = np.zeros((population_size, population_size))
    else:
        Meff_0_standardized = []
        Meff_1_standardized = []
        Meff_infinity_standardized = []
        jaccard_index = []
        is_same_day = []
        is_same_individual = []
        is_same_location = []

    for combination in list(combinations(range(population_size), 2)):
        x = combination[0]
        x_individual = 0
        x_individual_index = 0
        y = combination[1]
        y_individual = 0
        y_individual_index = 0
        for individual in range(len(individuals)):
            if x in individuals[individual]:
                x_individual = ds.individuals[individual]
                x_individual_index = x - individuals[individual][0]
            if y in individuals[individual]:
                y_individual = ds.individuals[individual]
                y_individual_index = y - individuals[individual][0]
        if not is_across_ds:
            sheet = x_individual.sheet
            is_combined_sheets = False
            split_sheets = False
        else:
            sheet = [x_individual.sheet, y_individual.sheet]
            is_combined_sheets = True
            split_sheets = True
        individual_data_temp = id.Individual(x_individual.wb, sheet, x_individual.name,
                                             [x_individual.x_dimension_starts[x_individual_index],
                                              y_individual.x_dimension_starts[y_individual_index]],
                                             [x_individual.x_dimension_finishes[x_individual_index],
                                              y_individual.x_dimension_finishes[y_individual_index]],
                                             [x_individual.y_dimension_starts[x_individual_index],
                                              y_individual.y_dimension_starts[y_individual_index]],
                                             [x_individual.y_dimension_finishes[x_individual_index],
                                              y_individual.y_dimension_finishes[y_individual_index]],
                                             [x_individual.data_points[x_individual_index],
                                              y_individual.data_points[y_individual_index]],
                                             is_combined_sheets=is_combined_sheets, split_sheets=split_sheets,
                                             is_simulation=is_simulation)
        if not is_pca:
            if not is_jaccard_index:
                Meff_0_standardized[x][y] = individual_data_temp.Meff_0_standardized
                Meff_0_standardized[y][x] = individual_data_temp.Meff_0_standardized
                Meff_1_standardized[x][y] = individual_data_temp.Meff_1_standardized
                Meff_1_standardized[y][x] = individual_data_temp.Meff_1_standardized
                Meff_infinity_standardized[x][y] = individual_data_temp.Meff_infinity_standardized
                Meff_infinity_standardized[y][x] = individual_data_temp.Meff_infinity_standardized
            else:
                jaccard_index[x][y] = individual_data_temp.jaccard_index
                jaccard_index[y][x] = individual_data_temp.jaccard_index
        else:
            Meff_0_standardized.append(individual_data_temp.Meff_0_standardized)
            Meff_1_standardized.append(individual_data_temp.Meff_1_standardized)
            Meff_infinity_standardized.append(individual_data_temp.Meff_infinity_standardized)
            jaccard_index.append(individual_data_temp.jaccard_index)
            if x_individual.sheet == y_individual.sheet:
                is_same_day.append(True)
                if x_individual.name == y_individual.name:
                    is_same_individual.append(True)
                else:
                    is_same_individual.append(False)
            else:
                is_same_day.append(False)
                is_same_individual.append(False)
            if x_individual.data_points[x_individual_index] == y_individual.data_points[y_individual_index]:
                is_same_location.append(True)
            else:
                is_same_location.append(False)
        if specific_lobe is not None and x_individual.data_points[x_individual_index] == specific_lobe and \
            y_individual.data_points[y_individual_index] == specific_lobe:
            Meff_0_values.append(individual_data_temp.Meff_0_standardized)
            Meff_1_values.append(individual_data_temp.Meff_1_standardized)
            Meff_infinity_values.append(individual_data_temp.Meff_infinity_standardized)
            jaccard_values.append(individual_data_temp.jaccard_index)

    if specific_lobe is not None:
        return Meff_0_values, Meff_1_values, Meff_infinity_values, jaccard_values

    if is_pca:
        return Meff_0_standardized, Meff_1_standardized, Meff_infinity_standardized, jaccard_index, is_same_day, \
               is_same_individual, is_same_location

    if not is_jaccard_index:

        df_Meff_0_standardized = pd.DataFrame(Meff_0_standardized, index=names,
                                              columns=names)
        df_Meff_1_standardized = pd.DataFrame(Meff_1_standardized, index=names,
                                              columns=names)
        df_Meff_infinity_standardized = pd.DataFrame(Meff_infinity_standardized, index=names,
                                                     columns=names)

        if solo is True:
            if ax is None:
                fig, ax = plt.subplots()
            mask = np.zeros_like(df_Meff_0_standardized)
            mask[np.triu_indices_from(mask)] = True

            if is_cube:
                if species_label == 'Pig Lung':
                    df_cube = df_Meff_0_standardized.iloc[37:45, 14:21]
                    fig, ax = plt.subplots(figsize=(7, 8))
                    sns.heatmap(df_cube, vmin=0, vmax=1, ax=ax, yticklabels=False, xticklabels=False, cmap='crest',
                                cbar=False)
                    plt.savefig(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 '
                                f'Results/Fig 4 Cube.tiff', dpi=300, format='tiff')
                    return
                else:
                    df_cube = df_Meff_0_standardized.iloc[26:28, 6:8]
                    fig, ax = plt.subplots(figsize=(2, 2))
                    sns.heatmap(df_cube, vmin=0, vmax=1, ax=ax, yticklabels=False, xticklabels=False, cmap='crest',
                                cbar=False)
                    plt.savefig(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 '
                                f'Results/Fig 2 Cube.tiff', dpi=300, format='tiff')
                    return
            if species_label == 'Pig Lung':
                ax.hlines(7, 0, 45, colors='w')
                ax.hlines(14, 0, 45, colors='w')
                ax.hlines(21, 0, 45, colors='w')
                ax.hlines(29, 0, 45, colors='w')
                ax.hlines(37, 0, 45, colors='w')

                ax.vlines(7, 0, 45, colors='w')
                ax.vlines(14, 0, 45, colors='w')
                ax.vlines(21, 0, 45, colors='w')
                ax.vlines(29, 0, 45, colors='w')
                ax.vlines(37, 0, 45, colors='w')
                # plt.plot([], [], ' ', label="1. L.Apical\n2. R.Apical\n3. L.Cardiac\n4. R.Cardiac\n5. L.Diaphr\n6. "
                #                             "R.Diaphr\n7. Intermed\n8. Nasal")
                # ax.legend(loc="upper right", handlelength=0, handletextpad=0)
            else:
                for animal in range(ds.number_of_individuals):
                    ax.hlines(animal*2, 0, 28, colors='w')
                    ax.vlines(animal * 2, 0, 28, colors='w')
                # plt.plot([], [], ' ', label="1. Lung\n2. Nasal Turbinate")
                # ax.legend(loc="upper right", handlelength=0, handletextpad=0)
            with sns.axes_style("white"):
                # df_Meff_0_standardized.to_excel(f'/Users/maxbagga/Desktop/OneDrive - Emory University/Box/Mallard Phylogeny Data/Pig Lung and Nasal Day 5 Results/Beta Diversity/Beta Diversity Matrix df.xlsx')
                ax = sns.heatmap(df_Meff_0_standardized, vmin=0, vmax=1, ax=ax, cbar_kws={'label': 'Beta Diversity',
                                                                                          'orientation':'horizontal'},
                                 yticklabels=False, xticklabels=False, mask=mask, cmap='crest')
                ax.figure.axes[-1].xaxis.label.set_size(24)
                for collection in ax.collections:
                    if collection is not None:
                        cbar = collection.colorbar
                cbar.ax.tick_params(labelsize=24)
            if species_label == 'Pig Lung':
                if excel_writer is not None:
                    df_Meff_0_standardized.to_excel(excel_writer, sheet_name='Fig 5A')
                plt.savefig(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 '
                            f'Results/Fig 5A {species_label} Beta Diversity Matrix.tiff', dpi=300, format='tiff')
            else:
                if excel_writer is not None:
                    df_Meff_0_standardized.to_excel(excel_writer, sheet_name='Fig 3A')
                pass
                # plt.savefig(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 '
                #             f'Results/Fig 2 {species_label} Beta Diversity Matrix.tiff', dpi=300, format='tiff')
        else:
            fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(30, 24 / 3))
            sns.heatmap(df_Meff_0_standardized, vmin=0, vmax=1, ax=ax1, cbar_kws={'label': 'Beta Diversity'})
            sns.heatmap(df_Meff_1_standardized, vmin=0, vmax=1, ax=ax2, cbar_kws={'label': 'Beta Diversity'})
            sns.heatmap(df_Meff_infinity_standardized, vmin=0, vmax=1, ax=ax3, cbar_kws={'label': 'Beta Diversity'})
            ax1.set_title('q = 0', fontsize=30)
            ax2.set_title('q = 1', fontsize=30)
            ax3.set_title('q = ∞', fontsize=30)
        # if not is_across_ds:
        #     if not is_simulation:
        #         plt.savefig(f'/Users/maxbagga/Desktop/OneDrive - Emory University/Box/Mallard Phylogeny Data/{ds.species_name} Results/Beta Diversity/'
        #                     f'{ds.species_name} Beta Diversity Matrix Across Individuals.png')
        #     else:
        #         plt.savefig(f'/Users/maxbagga/Desktop/OneDrive - Emory University/Box/Mallard Phylogeny Data/{ds.species_name} Results/Beta Diversity/'
        #                     f'{ds.species_name} Beta Diversity Matrix Across Individuals Simulated.png')
        # else:
        #     if not is_simulation:
        #         plt.savefig(f'/Users/maxbagga/Desktop/OneDrive - Emory University/Box/Mallard Phylogeny Data/{ds.species_name[0]} Results/Beta Diversity/'
        #                     f'{ds.species_name[0]} {ds.species_name[1]} Beta Diversity Matrix Across Individuals.png')
        #     else:
        #         plt.savefig(f'/Users/maxbagga/Desktop/OneDrive - Emory University/Box/Mallard Phylogeny Data/{ds.species_name[0]} Results/Beta Diversity/'
        #                     f'{ds.species_name[0]} {ds.species_name[1]} Beta Diversity Matrix Across Individuals '
        #                     f'Simulated.png')
    else:
        fig, ax1 = plt.subplots(ncols=1, figsize=(10, 24 / 3))

        df_jaccard_index = pd.DataFrame(jaccard_index, index=names, columns=names)

        sns.heatmap(df_jaccard_index, vmin=0, vmax=1, ax=ax1, cbar_kws={'label': 'Jaccard Index'})

        ax1.set_title('Jaccard Index', fontsize=30)

        if not is_across_ds:
            if not is_simulation:
                plt.savefig(f'/Users/maxbagga/Desktop/OneDrive - Emory University/Box/Mallard Phylogeny Data/{ds.species_name} Results/Jaccard Index/'
                            f'{ds.species_name} Jaccard Index Matrix Across Individuals.png')
            else:
                plt.savefig(f'/Users/maxbagga/Desktop/OneDrive - Emory University/Box/Mallard Phylogeny Data/{ds.species_name} Results/Jaccard Index/'
                            f'{ds.species_name} Jaccard Index Matrix Across Individuals Simulated.png')
        else:
            if not is_simulation:
                plt.savefig(f'/Users/maxbagga/Desktop/OneDrive - Emory University/Box/Mallard Phylogeny Data/{ds.species_name[0]} Results/Jaccard Index/'
                            f'{ds.species_name[0]} {ds.species_name[1]} Jaccard Index Matrix Across Individuals.png')
            else:
                plt.savefig(f'/Users/maxbagga/Desktop/OneDrive - Emory University/Box/Mallard Phylogeny Data/{ds.species_name[0]} Results/Jaccard Index/'
                            f'{ds.species_name[0]} {ds.species_name[1]} Jaccard Index Matrix Across Individuals '
                            f'Simulated.png')


def get_beta_diversity_matrix_across_ds(ds1, ds2, is_jaccard_index=False, is_simulation=False, specific_lobe=None,
                                        species_label=None, is_cube=False, excel_writer=None):
    """Creates heatmap matrix for the three beta diversities across datasets"""
    combined_ds = ds1 + ds2
    if specific_lobe is None:
        get_beta_diversity_matrix_across_indiviuals(combined_ds, True, is_jaccard_index=is_jaccard_index,
                                                    is_simulation=is_simulation, specific_lobe=None,
                                                    species_label=species_label, is_cube=is_cube,
                                                    excel_writer=excel_writer)
    else:
        return get_beta_diversity_matrix_across_indiviuals(combined_ds, True, is_jaccard_index=is_jaccard_index,
                                                           is_simulation=is_simulation, specific_lobe=specific_lobe,
                                                           is_cube=is_cube, excel_writer=excel_writer)


def get_t_test_across_individuals(ds1, ds2):
    Meff_0_values_1, Meff_1_values_1, Meff_infinity_values_1, names1 = get_beta_diversity_matrix(ds1, is_ttest=True)
    jaccard_values = get_beta_diversity_matrix(ds1, is_jaccard_index=True, is_ttest=True)
    Meff_0_values_2, Meff_1_values_2, Meff_infinity_values_2, names2 = get_beta_diversity_matrix(ds2, is_ttest=True)
    jaccard_values += get_beta_diversity_matrix(ds2, is_jaccard_index=True, is_ttest=True)
    Meff_0_values = Meff_0_values_1 + Meff_0_values_2
    Meff_1_values = Meff_1_values_1 + Meff_1_values_2
    Meff_infinity_values = Meff_infinity_values_1 + Meff_infinity_values_2
    names = names1 + names2

    size = len(Meff_0_values)
    Meff_0_standardized = np.zeros((size, size))
    Meff_1_standardized = np.zeros((size, size))
    Meff_infinity_standardized = np.zeros((size, size))
    jaccard_index = np.zeros((size, size))

    for combination in list(combinations(range(size), 2)):
        x = combination[0]
        y = combination[1]

        Meff_0_standardized[x][y] = ttest_ind(Meff_0_values[x], Meff_0_values[y])[1]
        Meff_0_standardized[y][x] = Meff_0_standardized[x][y]
        Meff_1_standardized[x][y] = ttest_ind(Meff_1_values[x], Meff_1_values[y])[1]
        Meff_1_standardized[y][x] = Meff_1_standardized[x][y]
        Meff_infinity_standardized[x][y] = ttest_ind(Meff_infinity_values[x], Meff_infinity_values[y])[1]
        Meff_infinity_standardized[y][x] = Meff_infinity_standardized[x][y]
        jaccard_index[x][y] = ttest_ind(jaccard_values[x], jaccard_values[y])[1]
        jaccard_index[y][x] = jaccard_index[x][y]

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, figsize=(40, 32 / 4))

    df_Meff_0_standardized = pd.DataFrame(Meff_0_standardized, index=names,columns=names)
    df_Meff_1_standardized = pd.DataFrame(Meff_1_standardized, index=names,columns=names)
    df_Meff_infinity_standardized = pd.DataFrame(Meff_infinity_standardized, index=names,columns=names)
    df_jaccard_index = pd.DataFrame(jaccard_index, index=names, columns=names)

    sns.heatmap(df_Meff_0_standardized, vmin=0, vmax=1, ax=ax1, cbar_kws={'label': 'p-Value'})
    sns.heatmap(df_Meff_1_standardized, vmin=0, vmax=1, ax=ax2, cbar_kws={'label': 'p-Value'})
    sns.heatmap(df_Meff_infinity_standardized, vmin=0, vmax=1, ax=ax3, cbar_kws={'label': 'p-Value'})
    sns.heatmap(df_jaccard_index, vmin=0, vmax=1, ax=ax4, cbar_kws={'label': 'p-Value'})

    ax1.set_title('q = 0', fontsize=30)
    ax2.set_title('q = 1', fontsize=30)
    ax3.set_title('q = ∞', fontsize=30)
    ax4.set_title('Jaccard Index', fontsize=30)

    plt.savefig(f'/Users/maxbagga/Desktop/Box/Mallard Phylogeny Data/{ds1.species_name} Results/Beta Diversity/'
                f'{ds1.species_name} t-Test Matrix Across Individuals.png')


def get_t_test_across_lobes(ds1, ds2):
    Meff_0_values = []
    Meff_1_values = []
    Meff_infinity_values = []
    jaccard_values = []
    for data_point in ds1.data_points[0]:
        print(data_point)
        Meff_0_values_1, Meff_1_values_1, Meff_infinity_values_1, jaccard_values_1 = \
            get_beta_diversity_matrix_across_ds(ds1, ds2, is_jaccard_index=False, specific_lobe=data_point)
        Meff_0_values.append(Meff_0_values_1)
        Meff_1_values.append(Meff_1_values_1)
        Meff_infinity_values.append(Meff_infinity_values)
        jaccard_values.append(jaccard_values_1)

    names = ds1.data_points[0]
    size = len(ds1.data_points[0])
    Meff_0_standardized = np.zeros((size, size))
    Meff_1_standardized = np.zeros((size, size))
    Meff_infinity_standardized = np.zeros((size, size))
    jaccard_index = np.zeros((size, size))

    for combination in list(combinations(range(size), 2)):
        x = combination[0]
        y = combination[1]

        Meff_0_standardized[x][y] = ttest_ind(Meff_0_values[x], Meff_0_values[y])[1]
        Meff_0_standardized[y][x] = Meff_0_standardized[x][y]
        Meff_1_standardized[x][y] = ttest_ind(Meff_1_values[x], Meff_1_values[y])[1]
        Meff_1_standardized[y][x] = Meff_1_standardized[x][y]
        Meff_infinity_standardized[x][y] = ttest_ind(Meff_infinity_values[x], Meff_infinity_values[y])[1]
        Meff_infinity_standardized[y][x] = Meff_infinity_standardized[x][y]
        jaccard_index[x][y] = ttest_ind(jaccard_values[x], jaccard_values[y])[1]
        jaccard_index[y][x] = jaccard_index[x][y]

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, figsize=(40, 32 / 4))

    df_Meff_0_standardized = pd.DataFrame(Meff_0_standardized, index=names, columns=names)
    df_Meff_1_standardized = pd.DataFrame(Meff_1_standardized, index=names, columns=names)
    df_Meff_infinity_standardized = pd.DataFrame(Meff_infinity_standardized, index=names, columns=names)
    df_jaccard_index = pd.DataFrame(jaccard_index, index=names, columns=names)

    sns.heatmap(df_Meff_0_standardized, vmin=0, vmax=1, ax=ax1, cbar_kws={'label': 'p-Value'})
    sns.heatmap(df_Meff_1_standardized, vmin=0, vmax=1, ax=ax2, cbar_kws={'label': 'p-Value'})
    sns.heatmap(df_Meff_infinity_standardized, vmin=0, vmax=1, ax=ax3, cbar_kws={'label': 'p-Value'})
    sns.heatmap(df_jaccard_index, vmin=0, vmax=1, ax=ax4, cbar_kws={'label': 'p-Value'})

    ax1.set_title('q = 0', fontsize=30)
    ax2.set_title('q = 1', fontsize=30)
    ax3.set_title('q = ∞', fontsize=30)
    ax4.set_title('Jaccard Index', fontsize=30)

    plt.savefig(f'/Users/maxbagga/Desktop/Box/Mallard Phylogeny Data/{ds1.species_name} Results/Beta Diversity/'
                f'{ds1.species_name} t-Test Matrix Across Lobes.png')


def PCA_three_factors(ds1, ds2, is_simulation=False):
    """PCA analysis on three factors: location, pig, and day on pig lung data"""
    combined_ds = ds1 + ds2
    Meff_0_standardized, Meff_1_standardized, Meff_infinity_standardized, jaccard_index, is_same_day, \
        is_same_individual, is_same_location = \
        get_beta_diversity_matrix_across_indiviuals(combined_ds, True, is_pca=True, is_simulation=is_simulation)
    df1 = pd.DataFrame({'Meff 0': Meff_0_standardized, 'Meff 1': Meff_1_standardized,
                        'Meff Infinity': Meff_infinity_standardized, 'Jaccard Index': jaccard_index})
    targets_df = pd.DataFrame({'Day': is_same_day, 'Individual': is_same_individual, 'Location': is_same_location})

    inputs = StandardScaler().fit_transform(df1)


    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(inputs)
    principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'])

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(30, 24/3))
    ax1.set_xlabel('Principal Component 1', fontsize=15)
    ax1.set_ylabel('Principal Component 2', fontsize=15)
    ax1.set_title('Same Day', fontsize=20)
    ax2.set_xlabel('Principal Component 1', fontsize=15)
    ax2.set_ylabel('Principal Component 2', fontsize=15)
    ax2.set_title('Same Individual', fontsize=20)
    ax3.set_xlabel('Principal Component 1', fontsize=15)
    ax3.set_ylabel('Principal Component 2', fontsize=15)
    ax3.set_title('Same Location', fontsize=20)

    targets = [True, False]
    colors = ['r', 'b']

    for target, color in zip(targets, colors):
        indicesToKeep = targets_df['Day'] == target
        ax1.scatter(principalDf.loc[indicesToKeep, 'principal component 1'],
                    principalDf.loc[indicesToKeep, 'principal component 2'], c=color, s=25)

        indicesToKeep = targets_df['Individual'] == target
        ax2.scatter(principalDf.loc[indicesToKeep, 'principal component 1'],
                    principalDf.loc[indicesToKeep, 'principal component 2'], c=color, s=25)

        indicesToKeep = targets_df['Location'] == target
        ax3.scatter(principalDf.loc[indicesToKeep, 'principal component 1'],
                    principalDf.loc[indicesToKeep, 'principal component 2'], c=color, s=25)

    if not is_simulation:
        plt.savefig('PCA Pig Lung Days 3 and 5')
    else:
        plt.savefig('PCA Pig Lung Days 3 and 5 Simulated')
    plt.show()
    print(pca.explained_variance_ratio_)


def reshuffle_analysis(ds, cycles=10, reshuffle_across_all=False, is_across_ds=None):
    if not is_across_ds:
        Meff_0_values, Meff_1_values, Meff_infinity_values, names = \
            get_beta_diversity_matrix(ds, is_jaccard_index=False, is_simulation=False, is_ttest=True)
        jaccard_values = get_beta_diversity_matrix(ds, is_jaccard_index=True, is_simulation=False, is_ttest=True)

        ds_richnesses = []
        ds_actual_lobe_richness = []
        ds_diversities = []
        ds_actual_lobe_diversity = []
        ds_beta_diversity_pigs_distributions = []
        ds_combination_names = []

        individual_richnesses = []
        individual_diversities = []

        count_Meff_0s = []
        count_Meff_1s = []
        count_Meff_infinitys = []
        count_jaccards = []

        avg_Meff_0s = []
        avg_Meff_1s = []
        avg_Meff_infinitys = []
        avg_jaccards = []

        actual_lobe_richness = [[day.genotype_analyzer.unique for day in individual.days] for individual in
                                ds.individuals]
        actual_lobe_diversity = [[day.genotype_analyzer.shannon_wiener for day in individual.days] for
                                 individual in ds.individuals]

        avg_richnesss = []
        avg_diversitys = []

        beta_diversity_pigs_percentiles = []
        beta_diversity_pigs_distributions = []
        combination_names_temp = []

        labels = [individual.name for individual in ds.individuals]

        number_of_occurences = 0

        for individual in range(len(ds.individuals)):
            number_of_occurences += 1
            count_Meff_0 = 0
            count_Meff_1 = 0
            count_Meff_infinity = 0
            count_jaccard = 0
            Meff_0 = Meff_0_values[individual]
            Meff_1 = Meff_1_values[individual]
            Meff_infinity = Meff_infinity_values[individual]
            jaccard = jaccard_values[individual]
            current_individual = ds.individuals[individual]

            avg_Meff_0 = []
            avg_Meff_1 = []
            avg_Meff_infinity = []
            avg_jaccard = []

            avg_richness = []
            avg_diversity = []

            q0_beta_diversities = []
            combination_names = []
            richnesses = [[] for x in range(ds.len_data_points[individual])]
            diversities = [[] for x in range(ds.len_data_points[individual])]

            for cycle in range(cycles):
                if cycle % 25 == 0:
                    print('Cycle:', cycle)

                days = current_individual.reshuffle_plaques()

                for day in range(len(days)):
                    if cycle == 0:
                        avg_richness.append(days[day].genotype_analyzer.unique)
                        avg_diversity.append(days[day].genotype_analyzer.shannon_wiener)
                    else:
                        avg_richness[day] += days[day].genotype_analyzer.unique
                        avg_diversity[day] += days[day].genotype_analyzer.shannon_wiener
                    richnesses[day].append(days[day].genotype_analyzer.unique)
                    diversities[day].append(days[day].genotype_analyzer.shannon_wiener)

                Meff_0_values_individual = []
                Meff_1_values_individual = []
                Meff_infinity_values_individual = []
                jaccard_values_individual = []

                combination_count = 0
                for combination in list(combinations(range(len(days)), 2)):
                    x = combination[0]
                    y = combination[1]
                    if reshuffle_across_all:
                        individual_data_temp = id.Individual(current_individual.wb, current_individual.sheet,
                                                             current_individual.name,
                                                             [current_individual.x_dimension_starts[x],
                                                              current_individual.x_dimension_starts[y]],
                                                             [current_individual.x_dimension_finishes[x],
                                                              current_individual.x_dimension_finishes[y]],
                                                             [current_individual.y_dimension_starts[x],
                                                              current_individual.y_dimension_starts[y]],
                                                             [current_individual.y_dimension_finishes[x],
                                                              current_individual.y_dimension_finishes[y]],
                                                             [current_individual.data_points[x],
                                                              current_individual.data_points[y]],
                                                             override_plaques=[days[x].plaques_string,
                                                                               days[y].plaques_string])
                    else:
                        individual_data_temp = id.Individual(current_individual.wb, current_individual.sheet,
                                                             current_individual.name,
                                                             [current_individual.x_dimension_starts[x],
                                                              current_individual.x_dimension_starts[y]],
                                                             [current_individual.x_dimension_finishes[x],
                                                              current_individual.x_dimension_finishes[y]],
                                                             [current_individual.y_dimension_starts[x],
                                                              current_individual.y_dimension_starts[y]],
                                                             [current_individual.y_dimension_finishes[x],
                                                              current_individual.y_dimension_finishes[y]],
                                                             [current_individual.data_points[x],
                                                              current_individual.data_points[y]])
                        individual_data_reshuffle = individual_data_temp.reshuffle_plaques()
                        individual_data_temp = id.Individual(current_individual.wb, current_individual.sheet,
                                                             current_individual.name,
                                                             [current_individual.x_dimension_starts[x],
                                                              current_individual.x_dimension_starts[y]],
                                                             [current_individual.x_dimension_finishes[x],
                                                              current_individual.x_dimension_finishes[y]],
                                                             [current_individual.y_dimension_starts[x],
                                                              current_individual.y_dimension_starts[y]],
                                                             [current_individual.y_dimension_finishes[x],
                                                              current_individual.y_dimension_finishes[y]],
                                                             [current_individual.data_points[x],
                                                              current_individual.data_points[y]],
                                                             override_plaques=[individual_data_reshuffle[0].plaques_string,
                                                                               individual_data_reshuffle[1].plaques_string])
                    Meff_0_values_individual.append(individual_data_temp.Meff_0_standardized)
                    Meff_1_values_individual.append(individual_data_temp.Meff_1_standardized)
                    Meff_infinity_values_individual.append(individual_data_temp.Meff_infinity_standardized)
                    jaccard_values_individual.append(individual_data_temp.jaccard_index)
                    if cycle == 0:
                        avg_Meff_0.append(individual_data_temp.Meff_0_standardized)
                        avg_Meff_1.append(individual_data_temp.Meff_1_standardized)
                        avg_Meff_infinity.append(individual_data_temp.Meff_infinity_standardized)
                        avg_jaccard.append(individual_data_temp.jaccard_index)
                        q0_beta_diversities.append([individual_data_temp.Meff_0_standardized])
                        combination_names.append(f'{current_individual.data_points[x]} '
                                                 f'{current_individual.data_points[y]}')
                    else:
                        avg_Meff_0[combination_count] += individual_data_temp.Meff_0_standardized
                        avg_Meff_1[combination_count] += individual_data_temp.Meff_1_standardized
                        avg_Meff_infinity[combination_count] += individual_data_temp.Meff_infinity_standardized
                        avg_jaccard[combination_count] += individual_data_temp.jaccard_index
                        q0_beta_diversities[combination_count].append(individual_data_temp.Meff_0_standardized)
                    combination_count += 1

                if ttest_ind(Meff_0, Meff_0_values_individual)[1] < 0.05:
                    count_Meff_0 += 1
                if ttest_ind(Meff_1, Meff_1_values_individual)[1] < 0.05:
                    count_Meff_1 += 1
                if ttest_ind(Meff_infinity, Meff_infinity_values_individual)[1] < 0.05:
                    count_Meff_infinity += 1
                if ttest_ind(jaccard, jaccard_values_individual)[1] < 0.05:
                    count_jaccard += 1

            count_Meff_0s.append(count_Meff_0)
            count_Meff_1s.append(count_Meff_1)
            count_Meff_infinitys.append(count_Meff_infinity)
            count_jaccards.append(count_jaccard)

            avg_Meff_0s.append([x / cycles for x in avg_Meff_0])
            avg_Meff_1s.append([x / cycles for x in avg_Meff_1])
            avg_Meff_infinitys.append([x / cycles for x in avg_Meff_infinity])
            avg_jaccards.append([x / cycles for x in avg_jaccard])

            avg_richnesss.append([x / cycles for x in avg_richness])
            avg_diversitys.append([x / cycles for x in avg_diversity])

            beta_diversity_percentiles = []
            beta_diversity_distributions = []
            for combination in range(len(combination_names)):
                beta_diversity_percentiles.append(percentileofscore(q0_beta_diversities[combination],
                                                                    Meff_0_values[individual][combination]))
                beta_diversity_distributions.append([q0_beta_diversities[combination],
                                                     Meff_0_values[individual][combination]])
            beta_diversity_pigs_percentiles.append(beta_diversity_percentiles)
            beta_diversity_pigs_distributions.append(beta_diversity_distributions)
            individual_richnesses.append(richnesses)
            individual_diversities.append(diversities)
            individual_combination_names.append(combination_names)
        ds_richnesses.append(individual_richnesses)
        ds_actual_lobe_richness.append(actual_lobe_richness)
        ds_diversities.append(individual_diversities)
        ds_actual_lobe_diversity.append(actual_lobe_diversity)
        ds_beta_diversity_pigs_distributions.append(beta_diversity_pigs_distributions)
        ds_combination_names.append(individual_combination_names)

        combination_names_temp = combination_names
    else:
        ds_richnesses = []
        ds_actual_lobe_richness = []
        ds_diversities = []
        ds_actual_lobe_diversity = []
        ds_beta_diversity_pigs_distributions = []
        ds_combination_names = []
        for ds in is_across_ds:
            Meff_0_values, Meff_1_values, Meff_infinity_values, names = \
                get_beta_diversity_matrix(ds, is_jaccard_index=False, is_simulation=False, is_ttest=True)
            jaccard_values = get_beta_diversity_matrix(ds, is_jaccard_index=True, is_simulation=False, is_ttest=True)

            count_Meff_0s = []
            count_Meff_1s = []
            count_Meff_infinitys = []
            count_jaccards = []

            avg_Meff_0s = []
            avg_Meff_1s = []
            avg_Meff_infinitys = []
            avg_jaccards = []

            actual_lobe_richness = [[day.genotype_analyzer.unique for day in individual.days] for individual in
                                    ds.individuals]
            actual_lobe_diversity = [[day.genotype_analyzer.shannon_wiener for day in individual.days] for
                                     individual in ds.individuals]

            individual_richnesses = []
            individual_diversities = []

            avg_richnesss = []
            avg_diversitys = []

            beta_diversity_pigs_percentiles = []
            beta_diversity_pigs_distributions = []
            individual_combination_names = []

            labels = [individual.name for individual in ds.individuals]

            number_of_occurences = 0
            for individual in range(len(ds.individuals)):
                number_of_occurences += 1
                count_Meff_0 = 0
                count_Meff_1 = 0
                count_Meff_infinity = 0
                count_jaccard = 0
                Meff_0 = Meff_0_values[individual]
                Meff_1 = Meff_1_values[individual]
                Meff_infinity = Meff_infinity_values[individual]
                jaccard = jaccard_values[individual]
                current_individual = ds.individuals[individual]

                avg_Meff_0 = []
                avg_Meff_1 = []
                avg_Meff_infinity = []
                avg_jaccard = []

                avg_richness = []
                avg_diversity = []

                q0_beta_diversities = []
                combination_names = []
                richnesses = [[] for x in range(ds.len_data_points[individual])]
                diversities = [[] for x in range(ds.len_data_points[individual])]

                for cycle in range(cycles):
                    if cycle % 25 == 0:
                        print('Cycle:', cycle)

                    days = current_individual.reshuffle_plaques()

                    for day in range(len(days)):
                        if cycle == 0:
                            avg_richness.append(days[day].genotype_analyzer.unique)
                            avg_diversity.append(days[day].genotype_analyzer.shannon_wiener)
                        else:
                            avg_richness[day] += days[day].genotype_analyzer.unique
                            avg_diversity[day] += days[day].genotype_analyzer.shannon_wiener
                        richnesses[day].append(days[day].genotype_analyzer.unique)
                        diversities[day].append(days[day].genotype_analyzer.shannon_wiener)

                    Meff_0_values_individual = []
                    Meff_1_values_individual = []
                    Meff_infinity_values_individual = []
                    jaccard_values_individual = []

                    combination_count = 0
                    for combination in list(combinations(range(len(days)), 2)):
                        x = combination[0]
                        y = combination[1]
                        if reshuffle_across_all:
                            individual_data_temp = id.Individual(current_individual.wb, current_individual.sheet,
                                                                 current_individual.name,
                                                                 [current_individual.x_dimension_starts[x],
                                                                  current_individual.x_dimension_starts[y]],
                                                                 [current_individual.x_dimension_finishes[x],
                                                                  current_individual.x_dimension_finishes[y]],
                                                                 [current_individual.y_dimension_starts[x],
                                                                  current_individual.y_dimension_starts[y]],
                                                                 [current_individual.y_dimension_finishes[x],
                                                                  current_individual.y_dimension_finishes[y]],
                                                                 [current_individual.data_points[x],
                                                                  current_individual.data_points[y]],
                                                                 override_plaques=[days[x].plaques_string,
                                                                                   days[y].plaques_string])
                        else:
                            individual_data_temp = id.Individual(current_individual.wb, current_individual.sheet,
                                                                 current_individual.name,
                                                                 [current_individual.x_dimension_starts[x],
                                                                  current_individual.x_dimension_starts[y]],
                                                                 [current_individual.x_dimension_finishes[x],
                                                                  current_individual.x_dimension_finishes[y]],
                                                                 [current_individual.y_dimension_starts[x],
                                                                  current_individual.y_dimension_starts[y]],
                                                                 [current_individual.y_dimension_finishes[x],
                                                                  current_individual.y_dimension_finishes[y]],
                                                                 [current_individual.data_points[x],
                                                                  current_individual.data_points[y]])
                            individual_data_reshuffle = individual_data_temp.reshuffle_plaques()
                            individual_data_temp = id.Individual(current_individual.wb, current_individual.sheet,
                                                                 current_individual.name,
                                                                 [current_individual.x_dimension_starts[x],
                                                                  current_individual.x_dimension_starts[y]],
                                                                 [current_individual.x_dimension_finishes[x],
                                                                  current_individual.x_dimension_finishes[y]],
                                                                 [current_individual.y_dimension_starts[x],
                                                                  current_individual.y_dimension_starts[y]],
                                                                 [current_individual.y_dimension_finishes[x],
                                                                  current_individual.y_dimension_finishes[y]],
                                                                 [current_individual.data_points[x],
                                                                  current_individual.data_points[y]],
                                                                 override_plaques=[
                                                                     individual_data_reshuffle[0].plaques_string,
                                                                     individual_data_reshuffle[1].plaques_string])
                        Meff_0_values_individual.append(individual_data_temp.Meff_0_standardized)
                        Meff_1_values_individual.append(individual_data_temp.Meff_1_standardized)
                        Meff_infinity_values_individual.append(individual_data_temp.Meff_infinity_standardized)
                        jaccard_values_individual.append(individual_data_temp.jaccard_index)
                        if cycle == 0:
                            avg_Meff_0.append(individual_data_temp.Meff_0_standardized)
                            avg_Meff_1.append(individual_data_temp.Meff_1_standardized)
                            avg_Meff_infinity.append(individual_data_temp.Meff_infinity_standardized)
                            avg_jaccard.append(individual_data_temp.jaccard_index)
                            q0_beta_diversities.append([individual_data_temp.Meff_0_standardized])
                            combination_names.append(f'{current_individual.data_points[x]} '
                                                     f'{current_individual.data_points[y]}')
                        else:
                            avg_Meff_0[combination_count] += individual_data_temp.Meff_0_standardized
                            avg_Meff_1[combination_count] += individual_data_temp.Meff_1_standardized
                            avg_Meff_infinity[combination_count] += individual_data_temp.Meff_infinity_standardized
                            avg_jaccard[combination_count] += individual_data_temp.jaccard_index
                            q0_beta_diversities[combination_count].append(individual_data_temp.Meff_0_standardized)
                        combination_count += 1

                    if ttest_ind(Meff_0, Meff_0_values_individual)[1] < 0.05:
                        count_Meff_0 += 1
                    if ttest_ind(Meff_1, Meff_1_values_individual)[1] < 0.05:
                        count_Meff_1 += 1
                    if ttest_ind(Meff_infinity, Meff_infinity_values_individual)[1] < 0.05:
                        count_Meff_infinity += 1
                    if ttest_ind(jaccard, jaccard_values_individual)[1] < 0.05:
                        count_jaccard += 1

                count_Meff_0s.append(count_Meff_0)
                count_Meff_1s.append(count_Meff_1)
                count_Meff_infinitys.append(count_Meff_infinity)
                count_jaccards.append(count_jaccard)

                avg_Meff_0s.append([x / cycles for x in avg_Meff_0])
                avg_Meff_1s.append([x / cycles for x in avg_Meff_1])
                avg_Meff_infinitys.append([x / cycles for x in avg_Meff_infinity])
                avg_jaccards.append([x / cycles for x in avg_jaccard])

                avg_richnesss.append([x / cycles for x in avg_richness])
                avg_diversitys.append([x / cycles for x in avg_diversity])

                beta_diversity_percentiles = []
                beta_diversity_distributions = []
                for combination in range(len(combination_names)):
                    beta_diversity_percentiles.append(percentileofscore(q0_beta_diversities[combination],
                                                                        Meff_0_values[individual][combination]))
                    beta_diversity_distributions.append([q0_beta_diversities[combination],
                                                         Meff_0_values[individual][combination]])
                beta_diversity_pigs_percentiles.append(beta_diversity_percentiles)
                beta_diversity_pigs_distributions.append(beta_diversity_distributions)
                individual_richnesses.append(richnesses)
                individual_diversities.append(diversities)
                individual_combination_names.append(combination_names)
            ds_richnesses.append(individual_richnesses)
            ds_actual_lobe_richness.append(actual_lobe_richness)
            ds_diversities.append(individual_diversities)
            ds_actual_lobe_diversity.append(actual_lobe_diversity)
            ds_beta_diversity_pigs_distributions.append(beta_diversity_pigs_distributions)
            ds_combination_names.append(individual_combination_names)

        if not reshuffle_across_all:
            fig, axs = plt.subplots(nrows=4, ncols=7, figsize=(70, 40))
            ds_combination_names_temp = []
            for ds_temp in ds_combination_names:
                ds_combination_names_temp += ds_temp
            combination_names_reference = max(ds_combination_names_temp, key=len)
            combination_names_reference_sorted = []
            for combo in combination_names_reference:
                words = [word.lower() for word in combo.split()]
                words.sort()
                combination_names_reference_sorted.append(words)
            colors = {'0': 'tab:brown',
                      '1': 'tab:olive',
                      '2': 'tab:gray'}
            is_first = True
            individual_count = 0
            values = [[] for n in range(len(combination_names_reference))]
            actual_values = [[] for n in range(len(combination_names_reference))]
            names = [[] for n in range(len(combination_names_reference))]
            actual_values_names = [[] for n in range(len(combination_names_reference))]
            for ds_number in range(len(is_across_ds)):
                ds = is_across_ds[ds_number]
                for individual in range(len(ds.individuals)):
                    for combination in range(len(ds_combination_names[ds_number][individual])):
                        words = [word.lower() for word in
                                 ds_combination_names[ds_number][individual][combination].split()]
                        words.sort()
                        index = combination_names_reference_sorted.index(words)
                        names[index] += [ds.individual_names[individual]] * len(
                            ds_beta_diversity_pigs_distributions[ds_number][individual][combination][0])
                        values[index] += ds_beta_diversity_pigs_distributions[ds_number][individual][combination][0]
                        actual_values[index] += \
                            [ds_beta_diversity_pigs_distributions[ds_number][individual][combination][1]]
                        actual_values_names[index] += [ds.individual_names[individual]]
                        # axs[int(index // 7), index % 7].scatter([individual_count]*len(
                        #     ds_beta_diversity_pigs_distributions[ds_number][individual][combination][0]),
                        #     ds_beta_diversity_pigs_distributions[ds_number][individual][combination][0],
                        #                                          color=colors.get(f'{individual}'))
                        # axs[int(index // 7), index % 7].hlines(
                        #     ds_beta_diversity_pigs_distributions[ds_number][individual][combination][1], -.25, 5.25,
                        #     color=colors.get(f'{individual}'), linestyles='dashed')
                        if (len(ds_combination_names[ds_number][individual]) == len(combination_names_reference)) \
                                and is_first:
                            axs[int(combination // 7), combination % 7].set_ylabel('Beta Diversity', fontsize=24)
                            axs[int(combination // 7), combination % 7].set_ylim(-.2, 1.2)
                            axs[int(combination // 7), combination % 7].set_xlim(-0.25, 5.25)
                            axs[int(combination // 7), combination % 7].set_xticks([0, 1, 2, 3, 4, 5])
                            axs[int(combination // 7), combination % 7].set_yticks([0, .2, .4, .6, .8, 1])
                            axs[int(combination // 7), combination % 7].tick_params(axis='y', labelsize=24)
                            individual_names = []
                            for ds_temp in is_across_ds:
                                individual_names += ds_temp.individual_names
                            axs[int(combination // 7), combination % 7].set_title(
                                combination_names_reference[combination], fontsize=24)
                    individual_count += 1
                    if (len(ds_combination_names[ds_number][individual]) == len(combination_names_reference)) and \
                            is_first: is_first = False
            pickle.dump([combination_names_reference, names, values, actual_values_names, actual_values],
                        open(f'{ds.species_name} Beta Diversity VPs {cycles}x1.p', "wb"))
            for combination in range(len(combination_names_reference)):
                sns.violinplot(orient='v', x=names[combination], y=values[combination], ax=axs[int(combination // 7),
                                                                                               combination % 7],
                               color='tab:gray')
                sns.swarmplot(orient='v', x=actual_values_names[combination][:], y=actual_values[combination],
                              ax=axs[int(combination // 7), combination % 7], color='r', size=20)
        else:
            lobe_names_reference = []
            max_richnesses_len = 0
            for ds_temp in range(len(ds_richnesses)):
                for individual in range(len(ds_richnesses[ds_temp])):
                    if len(ds_richnesses[ds_temp][individual]) > max_richnesses_len:
                        lobe_names_reference = is_across_ds[ds_temp].data_points[individual]
                        max_richnesses_len = len(ds_richnesses[ds_temp][individual])

            fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20, 10))

            richness_values = [[] for n in range(len(lobe_names_reference))]
            diversity_values = [[] for n in range(len(lobe_names_reference))]
            actual_richness_values = [[] for n in range(len(lobe_names_reference))]
            actual_diversity_values = [[] for n in range(len(lobe_names_reference))]
            names = [[] for n in range(len(lobe_names_reference))]
            actual_names = [[] for n in range(len(lobe_names_reference))]
            for ds_number in range(len(is_across_ds)):
                ds = is_across_ds[ds_number]
                for individual in range(len(ds.individuals)):
                    for lobe in range(len(ds.data_points[individual])):
                        words = ds.data_points[individual][lobe]
                        index = lobe_names_reference.index(words)
                        names[index] += [ds.individual_names[individual]] * \
                                        len(ds_richnesses[ds_number][individual][lobe])
                        richness_values[index] += ds_richnesses[ds_number][individual][lobe]
                        actual_richness_values[index] += [ds_actual_lobe_richness[ds_number][individual][lobe]]
                        actual_names[index] += [ds.individual_names[individual]]
                        # sns.violinplot(orient='v', data=richnesses[day], ax=axs[day][0])
                        # axs[day][0].hlines(actual_lobe_richness[individual][day], -1, 1, linestyles='dashed')
                        # axs[day][0].set_xlim(-1, 1)
                        # axs[day][0].set_ylim(0, 21)
                        # axs[day][0].set_xticks([0, 1])
                        # axs[day][0].set_yticks([0, 5, 10, 15, 20])
                        # axs[day][0].tick_params(axis='x', labelsize=0)

                        diversity_values[index] += ds_diversities[ds_number][individual][lobe]
                        actual_diversity_values[index] += [ds_actual_lobe_diversity[ds_number][individual][lobe]]
                        # sns.violinplot(orient='v', data=diversities[day], ax=axs[day][1])
                        # axs[day][1].hlines(actual_lobe_diversity[individual][day], -1, 1, linestyles='dashed')
                        # axs[day][1].set_xlim(-1, 1)
                        # axs[day][1].set_ylim(0, 4)
                        # axs[day][1].set_xticks([0, 1])
                        # axs[day][1].set_yticks([0, 1, 2, 3, 4])
                        # axs[day][1].tick_params(axis='x', labelsize=0)

            plot_richness = [item for sublist in richness_values for item in sublist]
            plot_diversity = [item for sublist in diversity_values for item in sublist]
            plot_names = [item for sublist in names for item in sublist]
            plot_richness_names = [[] for n in range(len(lobe_names_reference))]
            plot_diversity_names = [[] for n in range(len(lobe_names_reference))]
            plot_actual_richness = [[] for n in range(len(lobe_names_reference))]
            plot_actual_diversity = [[] for n in range(len(lobe_names_reference))]
            plot_actual_richness_names = [[] for n in range(len(lobe_names_reference))]
            plot_actual_diversity_names = [[] for n in range(len(lobe_names_reference))]
            plot_actual_individual_names = [[] for n in range(len(lobe_names_reference))]
            for lobe in range(len(lobe_names_reference)):
                plot_richness_names[lobe] += [lobe_names_reference[lobe]] * len(richness_values[lobe])
                plot_diversity_names[lobe] += [lobe_names_reference[lobe]] * len(diversity_values[lobe])

                plot_actual_richness[lobe] += actual_richness_values[lobe]
                plot_actual_diversity[lobe] += actual_diversity_values[lobe]
                plot_actual_richness_names[lobe] += [lobe_names_reference[lobe]] * len(actual_richness_values[lobe])
                plot_actual_diversity_names[lobe] += [lobe_names_reference[lobe]] * len(actual_diversity_values[lobe])
                plot_actual_individual_names[lobe] += actual_names[lobe]
            plot_richness_names2 = [item for sublist in plot_richness_names for item in sublist]
            plot_diversity_names2 = [item for sublist in plot_diversity_names for item in sublist]
            plot_actual_richness2 = [item for sublist in plot_actual_richness for item in sublist]
            plot_actual_diversity2 = [item for sublist in plot_actual_diversity for item in sublist]
            plot_actual_richness_names2 = [item for sublist in plot_actual_richness_names for item in sublist]
            plot_actual_diversity_names2 = [item for sublist in plot_actual_diversity_names for item in sublist]
            plot_actual_individual_names2 = [item for sublist in plot_actual_individual_names for item in sublist]

            pickle.dump([plot_names, plot_richness, plot_diversity, plot_actual_individual_names2,
                         plot_actual_richness2, plot_actual_richness_names2, plot_actual_diversity2,
                         plot_actual_diversity_names2],
                        open(f'{ds.species_name} Richness and Diversity VPs {cycles}x1.p', "wb"))

            sns.violinplot(orient='v', x=plot_names, y=plot_richness, ax=ax1, color='tab:gray')
            sns.violinplot(orient='v', x=plot_names, y=plot_diversity, ax=ax2, color='tab:gray')
            sns.swarmplot(orient='v', x=plot_actual_individual_names2, y=plot_actual_richness2, ax=ax1,
                          hue=plot_actual_richness_names2, size=10, palette="Paired")
            sns.swarmplot(orient='v', x=plot_actual_individual_names2, y=plot_actual_diversity2, ax=ax2,
                          hue=plot_actual_diversity_names2, size=10, palette="Paired")

            ax1.set_title('A', loc='left', fontsize=48, fontweight=500)
            ax2.set_title('B', loc='left', fontsize=48, fontweight=500)
            ax1.set_ylabel('Richness', fontsize=30)
            ax2.set_ylabel('Diversity', fontsize=30)
            ax2.get_legend().remove()

            # for ax, row in zip(axs[:, 0], lobe_names_reference):
            #     ax.set_ylabel(row, rotation=0, size='large')

            # plt.savefig('Trial 2.png')

    return
    averages = [avg_Meff_0s, avg_Meff_1s, avg_Meff_infinitys, avg_jaccards]

    x = np.arange(len(labels))
    width = 0.18  # the width of the bars

    fig, ax = plt.subplots(figsize=(10, 10))
    rects1 = ax.bar(x - 1.5*width, count_Meff_0s, width, label='q = 0', color='tab:brown')
    rects2 = ax.bar(x - .5*width, count_Meff_1s, width, label='q = 1', color='tab:olive')
    rects3 = ax.bar(x + .5*width, count_Meff_infinitys, width, label='q = ∞', color='tab:gray')
    rects4 = ax.bar(x + 1.5*width, count_jaccards, width, label='Jaccard', color='tab:cyan')

    ax.set_ylabel('Number of Significant Differences', fontsize=24)
    ax.set_xlabel('Individual', fontsize=24)
    ax.set_ylim(0, 1100)
    ax.set_xticks(x)
    ax.set_yticks([0, 200, 400, 600, 800, 1000])
    ax.tick_params(axis='y', labelsize=24)
    ax.set_xticklabels(labels, fontsize=24)
    ax.legend(fontsize=18)

    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=18)

    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)
    autolabel(rects4)

    pickle.dump([count_Meff_0s, count_Meff_1s, count_Meff_infinitys, count_jaccards],
                open(f'{ds.species_name} Reshuffled 1000x1.p', "wb"))
    plt.savefig(f'/Users/maxbagga/Desktop/Box/Mallard Phylogeny Data/{ds.species_name} Results/Beta Diversity/'
                f'{ds.species_name} Reshuffled {cycles}x1 Figure')

    colors = ['tab:brown', 'tab:olive', 'tab:gray', 'tab:cyan']

    for individual in range(len(labels)):
        fig, ax = plt.subplots(figsize=(10, 10))

        ax.set_ylabel('Diversity', fontsize=24)
        ax.set_ylim(-.05, 1.05)
        ax.set_xlim(-0.25, 3.25)
        ax.set_xticks(range(4))
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        ax.tick_params(axis='y', labelsize=24)
        ax.set_xticklabels(['q = 0', 'q = 1', 'q = ∞', 'Jaccard'], fontsize=24)

        for avg in range(len(averages)):
            xs = [avg]*len(averages[avg][individual])
            ax.scatter(xs, averages[avg][individual], color=colors[avg], s=120)

        plt.savefig(f'/Users/maxbagga/Desktop/Box/Mallard Phylogeny Data/{ds.species_name} Results/Beta Diversity/'
                    f'{ds.species_name} {labels[individual]} Reshuffled Average Beta Diversity {cycles}x1 Figure')

    pickle.dump(averages, open(f'{ds.species_name} Reshuffled Average Beta Diversity {cycles}x1.p', "wb"))

    colors = ['tab:brown', 'tab:olive', 'tab:gray']

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20, 10))
    for individual in range(len(actual_lobe_richness)):

        ax1.set_ylabel('Richness', fontsize=24)
        ax1.set_ylim(-.05, 21.05)
        ax1.set_xlim(-0.25, len(avg_richnesss[individual]) - 1 + 0.25)
        ax1.set_xticks(range(len(avg_richnesss[individual])))
        ax1.set_yticks([0, 5, 10, 15, 20])
        ax1.tick_params(axis='y', labelsize=24)
        ax1.set_xticklabels(ds.data_points[individual], fontsize=12)

        ax2.set_ylabel('Diversity', fontsize=24)
        ax2.set_ylim(-.15, 3.15)
        ax2.set_xlim(-0.25, len(avg_richnesss[individual]) - 1 + 0.25)
        ax2.set_xticks(range(len(avg_richnesss[individual])))
        ax2.set_yticks([0, .75, 1.5, 2.25, 3])
        ax2.tick_params(axis='y', labelsize=24)
        ax2.set_xticklabels(ds.data_points[individual], fontsize=12)

        for lobe in range(len(actual_lobe_richness[individual])):
            ax1.scatter(lobe, avg_richnesss[individual][lobe], color=colors[individual], s=120)
            ax2.scatter(lobe, avg_diversitys[individual][lobe], color=colors[individual], s=120)

            ax1.scatter(lobe, actual_lobe_richness[individual][lobe], color=colors[individual], s=120, marker='P')
            ax2.scatter(lobe, actual_lobe_diversity[individual][lobe], color=colors[individual], s=120, marker='P')

    plt.savefig(f'/Users/maxbagga/Desktop/Box/Mallard Phylogeny Data/{ds.species_name} Results/Beta Diversity/'
                f'{ds.species_name} Reshuffled Average Richness and Diversity {cycles}x1 Figure')
    plt.show()


