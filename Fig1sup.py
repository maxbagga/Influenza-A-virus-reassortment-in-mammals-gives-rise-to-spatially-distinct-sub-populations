import Fig2
import matplotlib.pyplot as plt
import pandas as pd


def fig1sup_develop(gp_ds, f_ds, p_ds):
    writer = pd.ExcelWriter(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 '
                            f'Results/Figure Data/Fig S1 Data.xlsx', engine='openpyxl')

    Fig2.fig2sup_develop(gp_ds, 5, 2, 'A', "Guinea Pig", excel_writer=writer)
    Fig2.fig2sup_develop(f_ds, 5, 2, 'B', "Ferret", excel_writer=writer)
    Fig2.fig2sup_develop(p_ds, 3, 2, 'C', "Pig", excel_writer=writer)
    writer.save()

