import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import BetaDiversityViolin as bdv
import pandas as pd


def get_richness_diversity_violin_plot(pickle_file, code, axs=None, excel_writer=None):
    plot_names, plot_richness, plot_diversity, plot_actual_individual_names2, plot_actual_richness2, \
    plot_actual_richness_names2, plot_actual_diversity2, plot_actual_diversity_names2 = pickle_file

    if excel_writer is not None:
        if code == 0:
            pd.DataFrame({'Lobe': plot_names, 'Richness': plot_richness, 'Diversity': plot_diversity}).to_excel(
                excel_writer
                ,
                sheet_name=
                'Fig 4G and'
                ' 4H Simula'
                'ted Values'
                )
            pd.DataFrame({'Lobe': plot_actual_richness_names2, 'Richness': plot_actual_richness2, 'Diversity':
                plot_actual_diversity2}).to_excel(excel_writer, sheet_name='Fig 4G and 4H Actual Values')
        else:
            pd.DataFrame({'Lobe': plot_names, 'Richness': plot_richness, 'Diversity': plot_diversity}).to_excel(
                excel_writer
                ,
                sheet_name=
                'Fig 2C and'
                ' 2D Simula'
                'ted Values'
                )
            pd.DataFrame({'Lobe': plot_actual_richness_names2, 'Richness': plot_actual_richness2, 'Diversity':
                plot_actual_diversity2}).to_excel(excel_writer, sheet_name='Fig 2C and 2D Actual Values')

    if code == 0:
        p79_index = plot_names.index('P79')
        p84_index = plot_names.index('P84')
        p73_index = plot_names.index('P73')

        plot_names = swapPositions(plot_names, 0, p79_index)
        plot_names = swapPositions(plot_names, 1, p84_index)
        plot_names = swapPositions(plot_names, 2, p73_index)

        plot_richness = swapPositions(plot_richness, 0, p79_index)
        plot_richness = swapPositions(plot_richness, 1, p84_index)
        plot_richness = swapPositions(plot_richness, 2, p73_index)

        plot_diversity = swapPositions(plot_diversity, 0, p79_index)
        plot_diversity = swapPositions(plot_diversity, 1, p84_index)
        plot_diversity = swapPositions(plot_diversity, 2, p73_index)

        p79_index = plot_actual_individual_names2.index('P79')
        p84_index = plot_actual_individual_names2.index('P84')
        p73_index = plot_actual_individual_names2.index('P73')

        plot_actual_individual_names2 = swapPositions(plot_actual_individual_names2, 0, p79_index)
        plot_actual_individual_names2 = swapPositions(plot_actual_individual_names2, 1, p84_index)
        plot_actual_individual_names2 = swapPositions(plot_actual_individual_names2, 2, p73_index)

        plot_actual_richness2 = swapPositions(plot_actual_richness2, 0, p79_index)
        plot_actual_richness2 = swapPositions(plot_actual_richness2, 1, p84_index)
        plot_actual_richness2 = swapPositions(plot_actual_richness2, 2, p73_index)

        plot_actual_diversity2 = swapPositions(plot_actual_diversity2, 0, p79_index)
        plot_actual_diversity2 = swapPositions(plot_actual_diversity2, 1, p84_index)
        plot_actual_diversity2 = swapPositions(plot_actual_diversity2, 2, p73_index)

        count = 0
        start = 0
        while(count < 6):
            ld_index = plot_actual_richness_names2.index('L.Diaphr')
            lc_index = plot_actual_richness_names2.index('L.Cardiac', start+1)

            plot_actual_richness_names2 = swapPositions(plot_actual_richness_names2, ld_index, lc_index)
            plot_actual_richness2 = swapPositions(plot_actual_richness2, ld_index, lc_index)
            plot_actual_diversity2 = swapPositions(plot_actual_diversity2, ld_index, lc_index)
            start = ld_index
            count += 1

        count = 0
        start = 0
        while(count < 6):
            rd_index = plot_actual_richness_names2.index('R.Diaphr')
            rc_index = plot_actual_richness_names2.index('R.Cardiac', start+1)

            plot_actual_richness_names2 = swapPositions(plot_actual_richness_names2, rd_index, rc_index)
            plot_actual_richness2 = swapPositions(plot_actual_richness2, rd_index, rc_index)
            plot_actual_diversity2 = swapPositions(plot_actual_diversity2, rd_index, rc_index)
            start = rd_index
            count += 1

        name_converter = {'P79': 'P1',
                          'P84': 'P2',
                          'P73': 'P3',
                          'P55': 'P4',
                          'P52': 'P5',
                          'P76': 'P6'}

        plot_names = [name_converter.get(name) for name in plot_names]
        plot_actual_individual_names2 = [name_converter.get(name) for name in plot_actual_individual_names2]

    if axs is None:
        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20, 10))
    else:
        ax1, ax2 = axs

    sns.violinplot(orient='v', x=plot_names, y=plot_richness, ax=ax1, color='tab:gray')
    sns.violinplot(orient='v', x=plot_names, y=plot_diversity, ax=ax2, color='tab:gray')

    percentile_booleans_richness, percentile_booleans_names_richness = \
        bdv.get_percentiles(plot_names, plot_richness, plot_actual_individual_names2, plot_actual_richness2,
                            starred=False)
    percentile_booleans_diversity, percentile_booleans_names_diversity = \
        bdv.get_percentiles(plot_names, plot_diversity, plot_actual_individual_names2, plot_actual_diversity2,
                            starred=False)

    x1 = plot_actual_individual_names2 #+ percentile_booleans_names_richness
    y1 = plot_actual_richness2 #+ percentile_booleans_richness
    z = plot_actual_richness_names2 #+ ['5th Percentile']*6
    x2 = plot_actual_individual_names2 #+ percentile_booleans_names_diversity
    y2 = plot_actual_diversity2 #+ percentile_booleans_diversity

    if code == 0:
        sns.swarmplot(orient='v', x=x1, y=y1, ax=ax1, hue=z, size=10, palette="Paired")
        sns.swarmplot(orient='v', x=x2, y=y2, ax=ax2, hue=z, size=10, palette="Paired")
    else:
        sns.swarmplot(orient='v', x=x1, y=y1, ax=ax1, hue=z, size=10, palette=sns.color_palette("tab10"))
        sns.swarmplot(orient='v', x=x2, y=y2, ax=ax2, hue=z, size=10, palette=sns.color_palette("tab10"))

    ticks = plt.gca().get_xticks()
    w = 0.3
    for actual_value in range(len(percentile_booleans_names_richness)):
        ax1.hlines(percentile_booleans_richness[actual_value], ticks[actual_value]-w, ticks[actual_value]+w, color='k')
    for actual_value in range(len(percentile_booleans_names_diversity)):
        ax2.hlines(percentile_booleans_diversity[actual_value], ticks[actual_value] - w, ticks[actual_value] + w,
                   color='k')

    percentile_booleans_richness, percentile_booleans_names_richness = \
        bdv.get_percentiles(plot_names, plot_richness, plot_actual_individual_names2, plot_actual_richness2,
                            starred=False, percentile_value=95)
    percentile_booleans_diversity, percentile_booleans_names_diversity = \
        bdv.get_percentiles(plot_names, plot_diversity, plot_actual_individual_names2, plot_actual_diversity2,
                            starred=False, percentile_value=95)

    for actual_value in range(len(percentile_booleans_names_richness)):
        ax1.hlines(percentile_booleans_richness[actual_value], ticks[actual_value]-w, ticks[actual_value]+w, color='k')
    for actual_value in range(len(percentile_booleans_names_diversity)):
        ax2.hlines(percentile_booleans_diversity[actual_value], ticks[actual_value] - w, ticks[actual_value] + w,
                   color='k')

    # sns.swarmplot(orient='v', x=percentile_booleans_names_richness, y=percentile_booleans_richness, ax=ax1, size=20, marker="*",alpha=.5)
    # sns.swarmplot(orient='v', x=percentile_booleans_names_diversity, y=percentile_booleans_diversity, ax=ax2, size=20, marker="*",alpha=.5)

    if code == 0:
        ax1.set_title('G', loc='left', fontsize=48, fontweight=500)
        ax2.set_title('H', loc='left', fontsize=48, fontweight=500)
        ax1.set_ylabel('Richness', fontsize=30)
        ax2.set_ylabel('Diversity', fontsize=30)
        ax1.set_ylim(-1, 21)
        ax1.set_yticks([0, 4, 8, 12, 16, 20])
        ax1.tick_params(axis='x', labelsize=24, rotation=45)
        ax1.tick_params(axis='y', labelsize=24)
        ax2.tick_params(axis='x', labelsize=24, rotation=45)
        ax2.tick_params(axis='y', labelsize=24)
        ax1.legend(fontsize=20, loc='upper left', ncol=3, handletextpad=0.05, columnspacing=0.1)
        ax2.get_legend().remove()
    else:
        ax1.set_title('C', loc='left', fontsize=48, fontweight=500)
        ax2.set_title('D', loc='left', fontsize=48, fontweight=500)
        ax1.set_ylabel('Richness', fontsize=30)
        ax2.set_ylabel('Diversity', fontsize=30)
        ax1.set_ylim(-1, 21)
        ax1.set_yticks([0, 4, 8, 12, 16, 20])
        ax1.tick_params(axis='x', labelsize=24, rotation=45)
        ax1.tick_params(axis='y', labelsize=24)
        ax2.tick_params(axis='x', labelsize=24, rotation=45)
        ax2.tick_params(axis='y', labelsize=24)
        ax1.legend(fontsize=20, loc='upper left', ncol=3, handletextpad=0.05, columnspacing=0.1)
        ax2.get_legend().remove()

    if code == 0:
        plt.savefig('/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/Fig 4'
                    'G and H.tiff', dpi=300, format='tiff')
    else:
        pass


def swapPositions(l, pos1, pos2):
    l[pos1], l[pos2] = l[pos2], l[pos1]
    return l
