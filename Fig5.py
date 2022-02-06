import RichnessDiversityViolin as rdv
import pickle
import BetaDiversityFigures as bdf
import BetaDiversityViolin as bdv
import pandas as pd


def fig5_develop(pig_lung_day3_ds, pig_lung_nasal_day5_ds):
    writer = pd.ExcelWriter(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 '
                            f'Results/Figure Data/Fig 5 Data.xlsx', engine='openpyxl')
    bdf.get_beta_diversity_matrix_across_ds(pig_lung_day3_ds, pig_lung_nasal_day5_ds, species_label='Pig Lung',
                                            excel_writer=writer)

    bdvp = pickle.load(open("Beta Diversity VPs 1000x1.p", "rb"))
    bdv.get_beta_diversity_violin_plot(bdvp, code=0, solo=True, excel_writer=writer)
    writer.save()
