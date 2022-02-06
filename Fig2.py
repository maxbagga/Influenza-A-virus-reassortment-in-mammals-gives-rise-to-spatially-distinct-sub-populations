import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib


def fig2_develop(ds, get_pf_swi=False, with_parentals=True, is_diversity=False, is_evenness=False, without_other=None,
                 get_remaining=False, additional_label=None, is_pie=False, is_bar=False, axs=None, target=None,
                 labels=None, code=2, is_fig6=False, is_fig2sup=False, tick_points=None, plot_points=None,
                 split_lines=False, is_x_label=False, excel_writer=False):
    """Develops figure 2 from previous study"""
    parent_frequencies = []
    uh_count = []
    swi = []

    parent_frequencies_df = []
    uh_count_df = []
    swi_df = []

    if get_remaining:
        remaining = []

    for animal in range(ds.number_of_individuals):
        gas = [day.genotype_analyzer for day in ds.individuals[animal].days]
        for ga1 in range(ds.individuals[animal].number_of_days):
            for ga2 in range(ds.individuals[animal].number_of_days):
                if ga1 != ga2:
                    gas[ga1], gas[ga2] = gas[ga1] + gas[ga2]

        pf_day_values = []
        uh_day_values = []
        swi_day_values = []

        if get_remaining:
            remaining_day_values = []

        for day in gas:
            H3N8 = 0
            H4N6 = 0
            for index in range(len(day.ids)):
                if day.ids[index] == "11111111":
                    H3N8 = day.appearances[index] / sum(day.appearances)
                elif day.ids[index] == "00000000":
                    H4N6 = day.appearances[index] / sum(day.appearances)
            pf_day_values.append(H3N8 + H4N6)
            if with_parentals:
                uh_day_values.append(day.unique)
            else:
                uh_day_values.append(day.unique_wop)
            if with_parentals and not is_evenness and without_other is None:
                swi_day_values.append(day.shannon_wiener)
            elif not with_parentals and not is_evenness and without_other is None:
                swi_day_values.append(day.shannon_wiener_wop)
            elif with_parentals and is_evenness and without_other is None:
                swi_day_values.append(day.evenness)
            elif not with_parentals and is_evenness and without_other is None:
                swi_day_values.append(day.evenness_wop)
            elif not with_parentals and not is_evenness and without_other is not None:
                swi_day_values.append(day.shannon_wiener_wop_woo)
            else:
                swi_day_values.append(day.evenness_wop_woo)
            if get_remaining and with_parentals and without_other is None:
                remaining_day_values.append(day.plaques)
            elif get_remaining and not with_parentals and without_other is None:
                remaining_day_values.append(day.plaques_wop)
            elif get_remaining and without_other is not None:
                remaining_day_values.append(day.plaques_wop_woo)
        parent_frequencies.append(pf_day_values)
        uh_count.append(uh_day_values)
        swi.append(swi_day_values)
        if len(pf_day_values) < ds.max_len_data_points:
            for difference in range(ds.max_len_data_points - len(pf_day_values)):
                pf_day_values.append(np.nan)
                uh_count.append(np.nan)
                swi.append(np.nan)
        parent_frequencies_df.append(pf_day_values)
        uh_count_df.append(uh_day_values)
        swi_df.append(swi_day_values)
        if get_remaining:
            remaining.append(remaining_day_values)

    if is_fig6:
        return parent_frequencies_df, uh_count_df, swi_df

    column_names = ds.individual_names
    df_pf = pd.DataFrame(parent_frequencies_df, columns=ds.individuals[ds.max_len_data_points_location].data_points).T
    df_uh = pd.DataFrame(uh_count_df, columns=ds.individuals[ds.max_len_data_points_location].data_points).T
    df_swi = pd.DataFrame(swi_df, columns=ds.individuals[ds.max_len_data_points_location].data_points).T
    df_pf.columns = column_names
    df_uh.columns = column_names
    df_swi.columns = column_names
    if get_remaining:
        df_remaining = pd.DataFrame(remaining, columns=ds.individuals[ds.max_len_data_points_location].data_points).T
        df_remaining.columns = column_names

    if get_pf_swi and not get_remaining:
        return df_pf, df_swi
    elif get_pf_swi and get_remaining:
        return df_pf, df_swi, df_remaining

    if labels is None:
        labels = ['A', 'B', 'C', 'D']

    stacks = []
    for animal in range(ds.number_of_individuals):
        if target is not None:
            if animal != target:
                continue
        if is_fig2sup:
            pass
        elif axs is None:
            if not ds.is_spatial:
                fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, figsize=(40, 20))
            else:
                fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(40, 20))
        else:
            ax1, ax2, ax3, ax4 = axs
        gas = [day.genotype_analyzer for day in ds.individuals[animal].days]
        for ga1 in range(ds.individuals[animal].number_of_days):
            for ga2 in range(ds.individuals[animal].number_of_days):
                if ga1 != ga2:
                    gas[ga1], gas[ga2] = gas[ga1] + gas[ga2]

        frequencies = []
        if with_parentals:
            for ga in gas:
                frequencies.append([app / sum(ga.appearances) for app in ga.appearances])
        else:
            for ga in gas:
                try:
                    frequencies.append([app / ga.plaques_wop for app in ga.appearances])
                except ZeroDivisionError:
                    frequencies.append([0 for app in ga.appearances])
        frequencies = pd.DataFrame(frequencies).T
        frequencies.columns = ds.individuals[animal].data_points

        H3N8_index = None
        H4N6_index = None
        for genotype in range(len(gas[0].ids)):
            if gas[0].ids[genotype] == "11111111":
                H3N8_index = genotype
            elif gas[0].ids[genotype] == "00000000":
                H4N6_index = genotype

        if H3N8_index is None and H4N6_index is None:
            pass
            style = 0
        elif H3N8_index is None:
            a, b = frequencies.iloc[0].copy(), frequencies.iloc[H4N6_index].copy()
            frequencies.iloc[0], frequencies.iloc[H4N6_index] = b, a
            if not with_parentals:
                frequencies.iloc[0] = [0] * ds.individuals[animal].number_of_days
            style = 1
        elif H4N6_index is None:
            a, b = frequencies.iloc[0].copy(), frequencies.iloc[H3N8_index].copy()
            frequencies.iloc[0], frequencies.iloc[H3N8_index] = b, a
            if not with_parentals:
                frequencies.iloc[0] = [0] * ds.individuals[animal].number_of_days
            style = 2
        else:
            a, b, c, d = frequencies.iloc[0].copy(), frequencies.iloc[1].copy(), frequencies.iloc[H3N8_index].copy(), \
                         frequencies.iloc[H4N6_index].copy()
            frequencies.iloc[0], frequencies.iloc[1], frequencies.iloc[H3N8_index], frequencies.iloc[H4N6_index] = c, \
                                                                                                                   d, \
                                                                                                                   a, b
            if not with_parentals:
                frequencies.iloc[0] = [0] * ds.individuals[animal].number_of_days
                frequencies.iloc[1] = [0] * ds.individuals[animal].number_of_days
            style = 3

        frequencies = frequencies.T

        if style == 0:
            frequencies.insert(0, '00000000', [0] * ds.individuals[animal].number_of_days, True)
            frequencies.insert(0, '11111111', [0] * ds.individuals[animal].number_of_days, True)
        elif style == 1:
            frequencies.insert(0, '11111111', [0] * ds.individuals[animal].number_of_days, True)
        elif style == 2:
            frequencies.insert(1, '00000000', [0] * ds.individuals[animal].number_of_days, True)

        if is_fig2sup:
            stacks.append(frequencies)
            continue
        if not is_pie and not is_bar:
            ax = frequencies.plot(kind='area', stacked=True, fontsize=30, legend=None,
                                  yticks=([0.2, 0.4, 0.6, 0.8, 1.0]), ax=ax1,
                                  xticks=([number for number in range(ds.individuals[animal].number_of_days)]))
            ax.set_ylabel('Genotype frequencies', fontsize=30)
            ax.set_xticklabels(ds.individuals[animal].data_points, fontsize=30)
            ax.margins(0)
            ax.set_ylim(0, 1)
            ax.set_title(labels[0], loc='left', fontsize=48, fontweight=500)
            if excel_writer is not None:
                if code == 0:
                    frequencies.to_excel(excel_writer, sheet_name='Fig 1A')
                if code == 1:
                    frequencies.to_excel(excel_writer, sheet_name='Fig 2A')
                if code == 2:
                    frequencies.to_excel(excel_writer, sheet_name='Fig 3A')
        elif is_pie:
            ax1.pie(ds.individuals[animal].appearances, shadow=True)
            ax1.set_title('A', loc='left', fontsize=48, fontweight=500)
            ax1.axis('equal')
        else:
            width = .35
            ax1.set_xlim(-0.5, 0.5 - width/2 + ds.individuals[ds.max_len_data_points_location].number_of_days - 1)
            ax1.set_xticks([number for number in range(ds.individuals[animal].number_of_days)])
            ax1.set_xticklabels(ds.individuals[animal].data_points, fontsize=30)
            ax1.set_title('A', loc='left', fontsize=48, fontweight=500)
            ax1.tick_params(axis='y', labelsize=30)
            ax1.set_ylabel('Number of genotypes', fontsize=30)
            for day in range(ds.individuals[animal].number_of_days):
                appearance_count = 0
                for appearance in ds.individuals[animal].days[day].appearances:
                    ax1.bar([day], [appearance], width, bottom=[appearance_count])
                    appearance_count += appearance

        if code == 0 or code == 1:
            styles = ['-']*6 + ['--']*4
        else:
            styles = ['-']*ds.number_of_individuals
        if tick_points is None:
            tick_points = ds.individuals[ds.max_len_data_points_location].number_of_days
            df_pf.index = np.arange(1, tick_points + 1)
            df_uh.index = np.arange(1, tick_points + 1)
            df_swi.index = np.arange(1, tick_points + 1)
        if plot_points is None:
            plot_points = range(ds.individuals[ds.max_len_data_points_location].number_of_days)

        if not split_lines:
            df_pf['Index'] = plot_points
            df_pf.set_index('Index', inplace=True)
            df_pf.index.name = ''
            df_uh['Index'] = plot_points
            df_uh.set_index('Index', inplace=True)
            df_uh.index.name = ''
            df_swi['Index'] = plot_points
            df_swi.set_index('Index', inplace=True)
            df_swi.index.name = ''
        else:
            df_pf1 = {}
            df_pf2 = {}
            df_uh1 = {}
            df_uh2 = {}
            df_swi1 = {}
            df_swi2 = {}
            count = 0
            for column in df_pf:
                if count < 6:
                    df_pf1[column] = df_pf[column].values.tolist()
                    df_uh1[column] = df_uh[column].values.tolist()
                    df_swi1[column] = df_swi[column].values.tolist()
                else:
                    df_pf2[column] = df_pf[column].values.tolist()
                    df_uh2[column] = df_uh[column].values.tolist()
                    df_swi2[column] = df_swi[column].values.tolist()
                count += 1
            df_pf_original = df_pf
            df_uh_original = df_uh
            df_swi_original = df_swi
            df_pf = pd.DataFrame(df_pf1)
            df_pf1 = pd.DataFrame(df_pf2)
            df_uh = pd.DataFrame(df_uh1)
            df_uh1 = pd.DataFrame(df_uh2)
            df_swi = pd.DataFrame(df_swi1)
            df_swi1 = pd.DataFrame(df_swi2)

            df_pf['Index'] = plot_points
            df_pf.set_index('Index', inplace=True)
            df_pf.index.name = ''

            df_pf1['Index'] = [1, 2, 3]
            df_pf1.set_index('Index', inplace=True)
            df_pf1.index.name = ''

            df_uh['Index'] = plot_points
            df_uh.set_index('Index', inplace=True)
            df_uh.index.name = ''

            df_uh1['Index'] = [1, 2, 3]
            df_uh1.set_index('Index', inplace=True)
            df_uh1.index.name = ''

            df_swi['Index'] = plot_points
            df_swi.set_index('Index', inplace=True)
            df_swi.index.name = ''

            df_swi1['Index'] = [1, 2, 3]
            df_swi1.set_index('Index', inplace=True)
            df_swi1.index.name = ''

        if split_lines:
            ax = df_pf.plot.line(marker='o', markersize=12, linewidth=4, fontsize=30, legend=None,
                                 xlim=(0.5, 0.5 + tick_points), xticks=([number + 1 for number in range(tick_points)]),
                                 ax=ax2, style=['-']*6)
            df_pf1.plot.line(marker='o', markersize=12, linewidth=4, fontsize=30, legend=None,
                             xlim=(0.5, 0.5 + tick_points), xticks=([number + 1 for number in range(tick_points)]),
                             ax=ax2, style=['--'] * 4)
        else:
            ax = df_pf.plot.line(marker='o', markersize=12, linewidth=4, fontsize=30, legend=None,
                                 xlim=(0.5, 0.5 + tick_points), xticks=([number + 1 for number in range(tick_points)]),
                                 ax=ax2, style=styles)
        ax.set_ylabel('Parental genotype frequencies', fontsize=30)
        if is_x_label:
            ax.set_xlabel('Days post-inoculation', fontsize=30)
        if tick_points is None:
            ax.set_xticklabels(ds.individuals[ds.max_len_data_points_location].data_points, fontsize=30)
        else:
            ax.set_xticklabels(range(1, 6), fontsize=30)
            # ax.set_xticklabels([f'Day {a+1}' for a in range(tick_points)], fontsize=30)
        ax.margins(0)
        ax.set_ylim(-.05, 1.05)
        ax.set_title(labels[1], loc='left', fontsize=48, fontweight=500)

        if split_lines:
            ax = df_uh.plot.line(marker='o', markersize=12, linewidth=4, fontsize=30, legend=None,
                                 xlim=(0.5, 0.5 + tick_points), ylim=(0, 22), yticks=([0, 5, 10, 15, 20]), ax=ax3,
                                 xticks=([number + 1 for number in range(tick_points)]), style=['-']*6)
            df_uh1.plot.line(marker='o', markersize=12, linewidth=4, fontsize=30, legend=None,
                             xlim=(0.5, 0.5 + tick_points), ylim=(0, 22), yticks=([0, 5, 10, 15, 20]), ax=ax3,
                             xticks=([number + 1 for number in range(tick_points)]), style=['--'] * 4)
        else:
            ax = df_uh.plot.line(marker='o', markersize=12, linewidth=4, fontsize=30, legend=None,
                                 xlim=(0.5, 0.5 + tick_points), ylim=(0, 22), yticks=([0, 5, 10, 15, 20]), ax=ax3,
                                 xticks=([number + 1 for number in range(tick_points)]), style=styles)
        if code == 0:
            ax.set_ylim(-1, 30)
            ax.legend(fontsize=20, ncol=3, loc='upper center', handlelength=3)
        if code == 1:
            ax.set_ylim(-1, 25)
            ax.legend(fontsize=20, ncol=3, loc='upper center', handlelength=3)
        if code == 2:
            ax.set_ylim(-1, 22)
            ax.legend(fontsize=20, ncol=3, loc='upper center', handlelength=3)
        ax.set_ylabel('Richness', fontsize=30)
        if is_x_label:
            ax.set_xlabel('Days post-inoculation', fontsize=30)
        if tick_points is None:
            ax.set_xticklabels(ds.individuals[ds.max_len_data_points_location].data_points, fontsize=30)
        else:
            ax.set_xticklabels(range(1, 6), fontsize=30)
            # ax.set_xticklabels([f'Day {a+1}' for a in range(tick_points)], fontsize=30)
        ax.margins(0)
        ax.set_title(labels[2], loc='left', fontsize=48, fontweight=500)

        if is_diversity:
            return df_swi

        if split_lines:
            ax = df_swi.plot.line(marker='o', markersize=12, linewidth=4, fontsize=30, legend=None,
                                  xlim=(0.5, 0.5 + tick_points), xticks=([number + 1 for number in range(tick_points)]),
                                  ax=ax4, style=['-']*6)
            df_swi1.plot.line(marker='o', markersize=12, linewidth=4, fontsize=30, legend=None,
                              xlim=(0.5, 0.5 + tick_points), xticks=([number + 1 for number in range(tick_points)]),
                              ax=ax4, style=['--'] * 4)
        else:
            ax = df_swi.plot.line(marker='o', markersize=12, linewidth=4, fontsize=30, legend=None,
                                  xlim=(0.5, 0.5 + tick_points), xticks=([number + 1 for number in range(tick_points)]),
                                  ax=ax4, style=styles)
        ax.set_ylabel('Diversity', fontsize=30)
        if is_x_label:
            ax.set_xlabel('Days post-inoculation', fontsize=30)
        if tick_points is None:
            ax.set_xticklabels(ds.individuals[ds.max_len_data_points_location].data_points, fontsize=30)
        else:
            ax.set_xticklabels(range(1, 6), fontsize=30)
            # ax.set_xticklabels([f'Day {a+1}' for a in range(tick_points)], fontsize=30)
        ax.margins(0)
        ax.set_ylim(-.15, 3.5)
        ax.set_title(labels[3], loc='left', fontsize=48, fontweight=500)

        if target is not None:
            if not split_lines:
                return df_pf, df_uh, df_swi
            else:
                return df_pf_original, df_pf, df_pf1, df_uh_original, df_uh, df_uh1, df_swi_original, df_swi, df_swi1,
    return stacks
        # if not with_parentals:
        #     if additional_label is None:
        #         plt.savefig(f'/Users/maxbagga/Box/Mallard Phylogeny Data/{ds.species_name} Results/Figure 2 Without '
        #                     f'Parentals/{ds.species_name} Figure 2 {ds.individuals[animal].name} Without Parentals.png')
        #     else:
        #         plt.savefig(f'/Users/maxbagga/Box/Mallard Phylogeny Data/{ds.species_name} Results/Figure 2 Without '
        #                     f'Parentals {additional_label}/{ds.species_name} Figure 2 {ds.individuals[animal].name} '
        #                     f'{additional_label} Without Parentals.png')
        # else:
        #     if additional_label is None:
        #         plt.savefig(f'/Users/maxbagga/Box/Mallard Phylogeny Data/{ds.species_name} Results/Figure 2/'
        #                     f'{ds.species_name} Figure 2 {ds.individuals[animal].name}.png')
        #     else:
        #         plt.savefig(f'/Users/maxbagga/Box/Mallard Phylogeny Data/{ds.species_name} Results/Figure 2 '
        #                     f'{additional_label}/{ds.species_name} Figure 2 {ds.individuals[animal].name} '
        #                     f'{additional_label}.png')


def fig2pie_charts_pig(ds, excel_writer=None):
    for animal in range(ds.number_of_individuals):
        fig, axs = plt.subplots(ncols=ds.individuals[animal].number_of_days, figsize=(
            ds.individuals[animal].number_of_days, 10))

        unique_genotypes_df = ds.individuals[animal].unique_genotypes_df.copy(deep=True)

        ugs_individual = unique_genotypes_df['Unique Genotypes'].tolist()

        H3N8_index = None
        H4N6_index = None
        for genotype in range(len(ugs_individual)):
            if ugs_individual[genotype] == "11111111":
                H3N8_index = genotype
            elif ugs_individual[genotype] == "00000000":
                H4N6_index = genotype

        colors_to_add = None
        if len(ugs_individual) <= 20:
            if H3N8_index is not None and H4N6_index is not None:
                temp = swapPositions(ugs_individual, 0, H3N8_index)
                swapped_ugs = swapPositions(temp, 2, H4N6_index)
            elif H3N8_index is not None:
                if len(ugs_individual) <= 19:
                    swapped_ugs = swapPositions(ugs_individual, 0, H3N8_index)
                    swapped_ugs.insert(2, '00000000')
                else:
                    swapped_ugs = swapPositions(ugs_individual, 0, H3N8_index)
                    swapped_ugs.insert(2, '00000000')
                    colors_to_add = 1
            elif H4N6_index is not None:
                if len(ugs_individual) <= 19:
                    swapped_ugs = swapPositions(ugs_individual, 1, H4N6_index)
                    swapped_ugs.insert(0, '11111111')
                else:
                    swapped_ugs = swapPositions(ugs_individual, 1, H4N6_index)
                    swapped_ugs.insert(0, '11111111')
                    colors_to_add = 1
            else:
                if len(ugs_individual) <= 18:
                    ugs_individual.insert(0, '11111111')
                    ugs_individual.insert(2, '00000000')
                    swapped_ugs = ugs_individual
                else:
                    ugs_individual.insert(0, '11111111')
                    ugs_individual.insert(2, '00000000')
                    swapped_ugs = ugs_individual
                    colors_to_add = 2
            if colors_to_add is None:
                cmap = matplotlib.cm.get_cmap('tab20')
                colors = cmap.colors
                new_colors = colors[:len(swapped_ugs)]
            else:
                cmap = matplotlib.cm.get_cmap('tab20')
                colors = cmap.colors
                new_colors = colors + matplotlib.cm.get_cmap('Pastel1').colors[:colors_to_add]
        else:
            if H3N8_index is not None and H4N6_index is not None:
                temp = swapPositions(ugs_individual, 20, H3N8_index)
                swapped_ugs = swapPositions(temp, 25, H4N6_index)
            elif H3N8_index is not None:
                swapped_ugs = swapPositions(ugs_individual, 20, H3N8_index)
                swapped_ugs.insert(25, '00000000')
            elif H4N6_index is not None:
                swapped_ugs = swapPositions(ugs_individual, 24, H4N6_index)
                swapped_ugs.insert(20, '11111111')
            else:
                ugs_individual.insert(20, '11111111')
                ugs_individual.insert(25, '00000000')
                swapped_ugs = ugs_individual
            new_colors = matplotlib.cm.get_cmap('tab20b').colors + \
                         matplotlib.cm.get_cmap('tab20c').colors[:len(swapped_ugs)-20]

        for lobe in range(ds.individuals[animal].number_of_days):
            df = ds.individuals[animal].days[lobe].unique_genotypes_df.copy(deep=True)
            ugs = df['Unique Genotypes'].tolist()
            appearances = df['Appearances'].tolist()
            plot_colors = []

            for ug in ugs:
                index = swapped_ugs.index(ug)
                plot_colors.append(new_colors[index])

            axs[lobe].pie(appearances, shadow=False, colors=plot_colors)
            axs[lobe].axis('equal')
            axs[lobe].set_title(ds.individuals[animal].days[lobe].day)

            if excel_writer is not None:
                excel_dic = {'P1': 'A', 'P2': 'B', 'P3': 'C', 'P4': 'D', 'P5': 'E', 'P6': 'F'}
                df['Appearances'].to_excel(excel_writer, sheet_name=f'Fig 4{excel_dic.get(ds.individuals[animal].name)}'
                                                                    f' {ds.individuals[animal].days[lobe].day}')


        plt.savefig(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/'
                    f'Pig Pie Charts/{ds.individuals[animal].name}.tiff', dpi=300, format='tiff')


def fig2pie_charts_ferret(ds, ncols, nrows, split_point=None, excel_writer=None):
    if split_point is None:
        fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(ncols * 5, 10))
        fig.suptitle('A', x='.05', fontsize=48, fontweight=500)
        lc = ds.number_of_individuals
        pie_charts(ds, lc, axs, False, excel_writer=excel_writer)
    else:
        lc = split_point
        lc1 = ds.number_of_individuals - split_point
        fig, axs = plt.subplots(ncols=lc//2, nrows=nrows*2, figsize=(ncols * 5, 20), constrained_layout=True)
        axs[0, 0].set_ylabel('Nasal', fontsize=30, rotation=0, labelpad=0)
        axs[1, 0].set_ylabel('Lung', fontsize=30, rotation=0, labelpad=0)
        axs[2, 0].set_ylabel('Nasal', fontsize=30, rotation=0, labelpad=0)
        axs[3, 0].set_ylabel('Lung', fontsize=30, rotation=0, labelpad=0)
        pie_charts(ds, lc, axs, split_point, 'A', 4, excel_writer=excel_writer)
        plt.clf()
        fig, axs = plt.subplots(ncols=lc1//2, nrows=nrows*2, figsize=(ncols * 5, 20), constrained_layout=True)
        axs[0, 0].set_ylabel('Nasal', fontsize=30, rotation=0, labelpad=0)
        axs[1, 0].set_ylabel('Lung', fontsize=30, rotation=0, labelpad=0)
        axs[2, 0].set_ylabel('Nasal', fontsize=30, rotation=0, labelpad=0)
        axs[3, 0].set_ylabel('Lung', fontsize=30, rotation=0, labelpad=0)
        pie_charts(ds, lc, axs, split_point, 'B', 3, ds.number_of_individuals, excel_writer=excel_writer)

    plt.savefig(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/Fig 2 '
                f'Ferret Pie Charts.tiff', dpi=300, format='tiff')


def pie_charts(ds, lc, axs, split_point=None, label=None, ncols=1, lc1=None, excel_writer=None):
    if lc1 is None:
        animals = range(lc)
    else:
        animals = range(lc, lc1)
    for animal in animals:
        unique_genotypes_df = ds.individuals[animal].unique_genotypes_df.copy(deep=True)

        ugs_individual = unique_genotypes_df['Unique Genotypes'].tolist()

        H3N8_index = None
        H4N6_index = None
        for genotype in range(len(ugs_individual)):
            if ugs_individual[genotype] == "11111111":
                H3N8_index = genotype
            elif ugs_individual[genotype] == "00000000":
                H4N6_index = genotype

        colors_to_add = None
        if len(ugs_individual) <= 20:
            if H3N8_index is not None and H4N6_index is not None:
                temp = swapPositions(ugs_individual, 0, H3N8_index)
                swapped_ugs = swapPositions(temp, 2, H4N6_index)
            elif H3N8_index is not None:
                if len(ugs_individual) <= 19:
                    swapped_ugs = swapPositions(ugs_individual, 0, H3N8_index)
                    swapped_ugs.insert(2, '00000000')
                else:
                    swapped_ugs = swapPositions(ugs_individual, 0, H3N8_index)
                    swapped_ugs.insert(2, '00000000')
                    colors_to_add = 1
            elif H4N6_index is not None:
                if len(ugs_individual) <= 19:
                    swapped_ugs = swapPositions(ugs_individual, 1, H4N6_index)
                    swapped_ugs.insert(0, '11111111')
                else:
                    swapped_ugs = swapPositions(ugs_individual, 1, H4N6_index)
                    swapped_ugs.insert(0, '11111111')
                    colors_to_add = 1
            else:
                if len(ugs_individual) <= 18:
                    ugs_individual.insert(0, '11111111')
                    ugs_individual.insert(2, '00000000')
                    swapped_ugs = ugs_individual
                else:
                    ugs_individual.insert(0, '11111111')
                    ugs_individual.insert(2, '00000000')
                    swapped_ugs = ugs_individual
                    colors_to_add = 2
            if colors_to_add is None:
                cmap = matplotlib.cm.get_cmap('tab20')
                colors = cmap.colors
                new_colors = colors[:len(swapped_ugs)]
            else:
                cmap = matplotlib.cm.get_cmap('tab20')
                colors = cmap.colors
                new_colors = colors + matplotlib.cm.get_cmap('Pastel1').colors[:colors_to_add]
        elif len(ugs_individual) < 26:
            if H3N8_index is not None and H4N6_index is not None:
                temp = swapPositions(ugs_individual, 0, H3N8_index)
                swapped_ugs = swapPositions(temp, 1, H4N6_index)
            elif H3N8_index is not None:
                swapped_ugs = swapPositions(ugs_individual, 0, H3N8_index)
                swapped_ugs.insert(1, '00000000')
            elif H4N6_index is not None:
                swapped_ugs = swapPositions(ugs_individual, 0, H4N6_index)
                swapped_ugs.insert(0, '11111111')
            else:
                ugs_individual.insert(0, '11111111')
                ugs_individual.insert(1, '00000000')
                swapped_ugs = ugs_individual
            if len(ugs_individual) >= 23:
                new_colors = tuple([matplotlib.cm.get_cmap('tab20c').colors[0]]) + \
                             tuple([matplotlib.cm.get_cmap('tab20c').colors[5]]) + \
                             matplotlib.cm.get_cmap('tab20b').colors + \
                             matplotlib.cm.get_cmap('tab20c').colors[20-len(swapped_ugs)+2:]
            else:
                new_colors = tuple([matplotlib.cm.get_cmap('tab20c').colors[0]]) + \
                             tuple([matplotlib.cm.get_cmap('tab20c').colors[5]]) + \
                             matplotlib.cm.get_cmap('tab20b').colors
        else:
            if H3N8_index is not None and H4N6_index is not None:
                temp = swapPositions(ugs_individual, 20, H3N8_index)
                swapped_ugs = swapPositions(temp, 25, H4N6_index)
            elif H3N8_index is not None:
                swapped_ugs = swapPositions(ugs_individual, 20, H3N8_index)
                swapped_ugs.insert(25, '00000000')
            elif H4N6_index is not None:
                swapped_ugs = swapPositions(ugs_individual, 24, H4N6_index)
                swapped_ugs.insert(20, '11111111')
            else:
                ugs_individual.insert(20, '11111111')
                ugs_individual.insert(25, '00000000')
                swapped_ugs = ugs_individual
            new_colors = matplotlib.cm.get_cmap('tab20b').colors + \
                         matplotlib.cm.get_cmap('tab20c').colors[:len(swapped_ugs)-20]
        if lc1 is not None:
            animal1 = animal - split_point
        else:
            animal1 = animal
        for lobe in range(ds.individuals[animal].number_of_days):
            df = ds.individuals[animal].days[lobe].unique_genotypes_df.copy(deep=True)
            ugs = df['Unique Genotypes'].tolist()
            appearances = df['Appearances'].tolist()
            plot_colors = []

            for ug in ugs:
                index = swapped_ugs.index(ug)
                plot_colors.append(new_colors[index])

            if split_point is not None:
                if animal1 // ncols < 1:
                    axs[1 - lobe, animal1].pie(appearances, shadow=False, colors=plot_colors)
                    axs[1 - lobe, animal1].axis('equal')
                    axs[0, animal1].set_title(ds.individuals[animal].name, fontsize=30, rotation=45)
                else:
                    axs[3 - lobe, animal1 % ncols].pie(appearances, shadow=False, colors=plot_colors)
                    axs[3 - lobe, animal1 % ncols].axis('equal')
                    axs[2, animal1 % ncols].set_title(ds.individuals[animal].name, fontsize=30, rotation=45)
                if excel_writer is not None:
                    df['Appearances'].to_excel(excel_writer, sheet_name=f'Fig 2{label} {ds.individuals[animal].name} '
                                                                        f'{ds.individuals[animal].days[lobe].day}')
            else:
                axs[1 - lobe, animal].pie(appearances, shadow=False, colors=plot_colors)
                axs[1 - lobe, animal].axis('equal')
                axs[0, animal].set_title(ds.individuals[animal].name, fontsize=30, rotation=45)
    if label is not None:
        plt.savefig(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/Fig '
                    f'2 Ferret Pie Charts {label}.tiff', dpi=300, format='tiff')
    else:
        plt.savefig(
            f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/Fig 2 '
            f'Ferret Pie Charts.tiff', dpi=300, format='tiff')


def fig2sup_develop(ds, cols, rows, suptitle, label, excel_writer=None):
    stacks = fig2_develop(ds, is_fig2sup=True)
    # stacks[6].to_excel('GP7f.xlsx')
    fig, axs = plt.subplots(ncols=cols, nrows=rows, figsize=(cols*10, rows*10), constrained_layout=True)
    count = 0
    for row in axs:
        for axtemp in row:
            ax = stacks[count].plot(kind='area', stacked=True, fontsize=30, legend=None,
                                    title=ds.individuals[count].name, yticks=([0.2, 0.4, 0.6, 0.8, 1.0]), ax=axtemp,
                                    xticks=([number for number in range(ds.individuals[count].number_of_days)]))
            ax.title.set_size(30)
            ax.set_ylabel('Genotype frequencies', fontsize=30)
            ax.set_xticklabels(ds.individuals[count].data_points, fontsize=30)
            ax.margins(0)
            ax.set_ylim(0, 1)
            # ax.set_title(labels[0], loc='left', fontsize=48, fontweight=500)
            # fig.suptitle(suptitle, x='.05', fontsize=48, fontweight=500)
            if excel_writer is not None:
                stacks[count].to_excel(excel_writer, sheet_name=f'Fig S1{suptitle} {ds.individuals[count].name}')
            count += 1

    plt.savefig(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/'
                f'Fig 1 supplementary {label}.tiff', dpi=300, format='tiff')


def swapPositions(l, pos1, pos2):
    l[pos1], l[pos2] = l[pos2], l[pos1]
    return l








