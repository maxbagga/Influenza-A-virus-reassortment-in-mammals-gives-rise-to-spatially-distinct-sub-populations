import BetaDiversityFigures as bdf
import pickle
import BetaDiversityViolin as bdv
import pandas as pd


def fig3sup_develop():
    writer = pd.ExcelWriter(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 '
                            f'Results/Figure Data/Fig S3 Data.xlsx', engine='openpyxl')

    bdvp = pickle.load(open("Beta Diversity VPs 1000x1.p", "rb"))
    bdv.get_beta_diversity_violin_plot(bdvp, code=0, excel_writer=writer)
    writer.save()
