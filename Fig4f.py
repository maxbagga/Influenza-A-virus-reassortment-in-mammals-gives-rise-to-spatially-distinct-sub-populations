import BetaDiversityFigures as bdf
import matplotlib.pyplot as plt
import pickle
import BetaDiversityViolin as bdv
import Fig2
import pandas as pd
import RichnessDiversityViolin as rdv


def fig4_develop(pig_lung_day3_ds, pig_lung_nasal_day5_ds):
    writer = pd.ExcelWriter(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 '
                            f'Results/Figure Data/Fig 4 Data.xlsx', engine='openpyxl')

    Fig2.fig2pie_charts_pig(pig_lung_day3_ds, excel_writer=writer)
    Fig2.fig2pie_charts_pig(pig_lung_nasal_day5_ds, excel_writer=writer)


    rdvp = pickle.load(open("Richness and Diversity VPs 1000x1.p", "rb"))
    rdv.get_richness_diversity_violin_plot(rdvp, 0, excel_writer=writer)
    writer.save()
