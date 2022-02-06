import BetaDiversityFigures as bdf
import pandas as pd


def fig2sup_develop(p3_ds, p5_ds, f_ds):

    bdf.get_beta_diversity_matrix(p3_ds, solo=True, species_label='Pig Lung Day 3', ncols=3, nrows=1, suptitle='A',
                                  code=1)
    bdf.get_beta_diversity_matrix(p5_ds, solo=True, species_label='Pig Lung and Nasal Day 5', ncols=3, nrows=1,
                                  suptitle='B', code=0)
    # bdf.get_beta_diversity_matrix(f_ds, solo=True, species_label='Ferret Lung and Nasal Turbinate', ncols=7, nrows=2,
    #                               suptitle='C', code=2)

