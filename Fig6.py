import Fig2
from scipy.stats import f_oneway
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def fig6_develop(gp_ds, f_ds, p_ds, axs=None, excel_writer=False):
    gp_ten_squared, gp_ten_fifth = gp_ds.get_split()
    gp_pf, gp_richness, gp_diversity = Fig2.fig2_develop(gp_ten_fifth, is_fig6=True)
    plot_gp_pf = [item for sublist in gp_pf for item in sublist]
    plot_gp_richness = [item for sublist in gp_richness for item in sublist]
    plot_gp_diversity = [item for sublist in gp_diversity for item in sublist]

    f_ten_squared, f_ten_fifth = f_ds.get_split()
    f_pf, f_richness, f_diversity = Fig2.fig2_develop(f_ten_fifth, is_fig6=True)
    plot_f_pf = [item for sublist in f_pf for item in sublist]
    plot_f_richness = [item for sublist in f_richness for item in sublist]
    plot_f_diversity = [item for sublist in f_diversity for item in sublist]

    p_pf, p_richness, p_diversity = Fig2.fig2_develop(p_ds, is_fig6=True)
    plot_p_pf = [item for sublist in p_pf for item in sublist]
    plot_p_richness = [item for sublist in p_richness for item in sublist]
    plot_p_diversity = [item for sublist in p_diversity for item in sublist]

    print('Parental Frequency ANOVA—')
    fvalue, pvalue = f_oneway(plot_gp_pf, plot_f_pf, plot_p_pf)
    print(f'All three: {pvalue}')
    fvalue, pvalue = f_oneway(plot_gp_pf, plot_f_pf)
    print(f'GP F: {pvalue}')
    fvalue, pvalue = f_oneway(plot_gp_pf,  plot_p_pf)
    print(f'GP P: {pvalue}')
    fvalue, pvalue = f_oneway(plot_f_pf, plot_p_pf)
    print(f'F P: {pvalue}')

    print('Richness ANOVA—')
    fvalue, pvalue = f_oneway(plot_gp_richness, plot_f_richness, plot_p_richness)
    print(f'All three: {pvalue}')
    fvalue, pvalue = f_oneway(plot_gp_richness, plot_f_richness)
    print(f'GP F: {pvalue}')
    fvalue, pvalue = f_oneway(plot_gp_richness,  plot_p_richness)
    print(f'GP P: {pvalue}')
    fvalue, pvalue = f_oneway(plot_f_richness, plot_p_richness)
    print(f'F P: {pvalue}')

    print('Diversity ANOVA—')
    fvalue, pvalue = f_oneway(plot_gp_diversity, plot_f_diversity, plot_p_diversity)
    print(f'All three: {pvalue}')
    fvalue, pvalue = f_oneway(plot_gp_diversity, plot_f_diversity)
    print(f'GP F: {pvalue}')
    fvalue, pvalue = f_oneway(plot_gp_diversity,  plot_p_diversity)
    print(f'GP P: {pvalue}')
    fvalue, pvalue = f_oneway(plot_f_diversity, plot_p_diversity)
    print(f'F P: {pvalue}')

    pfs = []
    richnesses = []
    diversities = []
    names = []
    for richness in range(len(plot_gp_richness)):
        pfs.append(plot_gp_pf[richness])
        richnesses.append(plot_gp_richness[richness])
        diversities.append(plot_gp_diversity[richness])
        names.append('GP')
    for richness in range(len(plot_f_richness)):
        pfs.append(plot_f_pf[richness])
        richnesses.append(plot_f_richness[richness])
        diversities.append(plot_f_diversity[richness])
        names.append('Ferret')
    for richness in range(len(plot_p_richness)):
        pfs.append(plot_p_pf[richness])
        richnesses.append(plot_p_richness[richness])
        diversities.append(plot_p_diversity[richness])
        names.append('Pig')

    if axs is None:
        fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(30, 10))
    else:
        ax1, ax2, ax3 = axs

    sns.violinplot(orient='v', x=names, y=pfs, ax=ax1, palette="tab10", cut=0)
    sns.violinplot(orient='v', x=names, y=richnesses, ax=ax2, palette="tab10", cut=0)
    sns.violinplot(orient='v', x=names, y=diversities, ax=ax3, palette="tab10", cut=0)

    if excel_writer:
        pd.DataFrame({'Animal': names, 'Parental Frequencies': pfs}).to_excel(excel_writer, sheet_name='Fig 1M')
        pd.DataFrame({'Animal': names, 'Richness': pfs}).to_excel(excel_writer, sheet_name='Fig 1N')
        pd.DataFrame({'Animal': names, 'Diversity': pfs}).to_excel(excel_writer, sheet_name='Fig 1O')

    ax1.set_title('M', loc='left', fontsize=48, fontweight=500)
    ax2.set_title('N', loc='left', fontsize=48, fontweight=500)
    ax3.set_title('O', loc='left', fontsize=48, fontweight=500)
    ax1.set_ylabel('Parental genotype frequencies', fontsize=30)
    ax2.set_ylabel('Richness', fontsize=30)
    ax3.set_ylabel('Diversity', fontsize=30)
    ax1.set_ylim(-.05, 1.05)
    ax2.set_ylim(-1, 22)
    ax3.set_ylim(-.15, 3.5)
    ax2.set_yticks([0, 5, 10, 15, 20])
    ax3.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
    ax1.tick_params(axis='x', labelsize=30)
    ax1.tick_params(axis='y', labelsize=30)
    ax2.tick_params(axis='x', labelsize=30)
    ax2.tick_params(axis='y', labelsize=30)
    ax3.tick_params(axis='x', labelsize=30)
    ax3.tick_params(axis='y', labelsize=30)

    plt.savefig('/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/'
                'Fig 6.tiff', dpi=300, format='tiff')

