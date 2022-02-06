import Fig2
import BetaDiversityFigures as bdf
import BetaDiversityViolin as bdv
from scipy.stats import ttest_rel
import pickle
import matplotlib.pyplot as plt
import pandas as pd
import RichnessDiversityViolin as rdv
import seaborn as sns


def fig2_develop(f_ds):
    writer = pd.ExcelWriter(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 '
                            f'Results/Figure Data/Fig 2 Data.xlsx', engine='openpyxl')
    Fig2.fig2pie_charts_ferret(f_ds, 14, 2, split_point=8, excel_writer=writer)
    # bdf.reshuffle_analysis(f_ds, cycles=1000, is_across_ds=[f_ds], reshuffle_across_all=True)

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(20, 20))

    parent_frequencies_df, f_richness, f_diversity = Fig2.fig2_develop(f_ds, is_fig6=True)

    lung_richness = []
    lung_diversity = []
    nasal_richness = []
    nasal_diversity = []

    for animal in range(len(f_richness)):
        lung_richness.append(f_richness[animal][0])
        nasal_richness.append(f_richness[animal][1])
        lung_diversity.append(f_diversity[animal][0])
        nasal_diversity.append(f_diversity[animal][1])

    print(ttest_rel(lung_richness, nasal_richness))
    print(ttest_rel(lung_diversity, nasal_diversity))

    richness = lung_richness + nasal_richness
    diversity = lung_diversity + nasal_diversity

    names = ['Lung'] * len(lung_richness) + ['Nasal Turbinate'] * len(nasal_richness)

    ax3 = axs[1, 0]
    ax4 = axs[1, 1]
    ax1 = axs[0, 0]
    ax2 = axs[0, 1]

    rdvp = pickle.load(open("Ferret Lung and Nasal Turbinate Richness and Diversity VPs 1000x1.p", "rb"))
    rdv.get_richness_diversity_violin_plot(rdvp, 1, axs=(ax3, ax4), excel_writer=writer)

    sns.violinplot(orient='v', x=names, y=richness, ax=ax1, palette="tab10")
    sns.violinplot(orient='v', x=names, y=diversity, ax=ax2, palette="tab10")

    pd.DataFrame({'Lobe': names, 'Richness': richness}).to_excel(writer, sheet_name='Fig 2E')
    pd.DataFrame({'Lobe': names, 'Diversity': diversity}).to_excel(writer, sheet_name='Fig 2F')

    writer.save()

    ax1.set_title('E', loc='left', fontsize=48, fontweight=500)
    ax2.set_title('F', loc='left', fontsize=48, fontweight=500)
    ax1.set_ylabel('Richness', fontsize=30)
    ax2.set_ylabel('Diversity', fontsize=30)
    ax1.set_ylim(-1, 21)
    ax1.set_yticks([0, 4, 8, 12, 16, 20])
    ax2.set_yticks([0.0, 1.0, 2.0, 3.0])
    ax1.tick_params(axis='x', labelsize=24)
    ax1.tick_params(axis='y', labelsize=24)
    ax2.tick_params(axis='x', labelsize=24)
    ax2.tick_params(axis='y', labelsize=24)
    ax1.plot([], [], ' ', label=f'p-value = {round(ttest_rel(lung_richness, nasal_richness)[1], 6)}')
    ax1.legend(loc="lower right", handlelength=0, handletextpad=0, fontsize=20)
    ax2.plot([], [], ' ', label=f'p-value = {round(ttest_rel(lung_diversity, nasal_diversity)[1], 4)}')
    ax2.legend(loc="lower right", handlelength=0, handletextpad=0, fontsize=20)

    pos1 = ax1.get_position()
    pos2 = ax2.get_position()
    pos3 = ax3.get_position()
    pos4 = ax4.get_position()
    ax1.set_position(pos3)
    ax2.set_position(pos4)
    ax3.set_position(pos1)
    ax4.set_position(pos2)

    plt.savefig('/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/Fig 2'
                '.tiff', dpi=300, format='tiff')
