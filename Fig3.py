import BetaDiversityFigures as bdf
import pickle
import RichnessDiversityViolin as rdv
import matplotlib.pyplot as plt
import Fig2
import seaborn as sns
from scipy.stats import ttest_rel
import BetaDiversityViolin as bdv
import pandas as pd


def fig3_develop(f_ds):
    writer = pd.ExcelWriter(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 '
                            f'Results/Figure Data/Fig 3 Data.xlsx', engine='openpyxl')
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(20, 10))

    bdf.get_beta_diversity_matrix_across_indiviuals(f_ds, species_label='Ferret', ax=ax1, excel_writer=writer)

    # bdf.reshuffle_analysis(f_ds, cycles=1000, is_across_ds=[f_ds], reshuffle_across_all=False)

    bdvp = pickle.load(open("Ferret Lung and Nasal Turbinate Beta Diversity VPs 1000x1.p", "rb"))
    bdv.get_beta_diversity_violin_plot(bdvp, 1, ax=ax2, excel_writer=writer)
    writer.save()

    ax1.set_title('A', loc='left', fontsize=48, fontweight=500)
    ax2.set_title('B', loc='left', fontsize=48, fontweight=500)

    plt.savefig('/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/'
                'Fig 3.tiff', dpi=300, format='tiff')

    bdf.get_beta_diversity_matrix_across_indiviuals(f_ds, species_label='Ferret', is_cube=True)



