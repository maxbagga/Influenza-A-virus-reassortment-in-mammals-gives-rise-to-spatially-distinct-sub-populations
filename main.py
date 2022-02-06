import SpeciesData as sd
import Fig2
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
import pickle
import random
import scipy
from scipy import stats
import seaborn as sns
import CreateSpecies as cs
import BetaDiversityFigures as bdf
import RichnessDiversityViolin as rdv
import BetaDiversityViolin as bdv
from collections import OrderedDict
from scipy.stats import f_oneway
import ANOVATests as at
import Fig1
import Fig2f
import Fig3
import Fig4f
import Fig5
import Fig6
import Fig1sup
import Fig2sup
import Fig3sup


guinea_pig_nasal_ds = cs.create_guinea_pig_nasal()
ferret_nasal_ds = cs.create_ferret_nasal()
pig_nasal_ds = cs.create_pig_nasal()
ferret_lung_nasal_ds = cs.create_ferret_lung_nasal()
pig_lung_day3_ds = cs.create_pig_lung_day3()
pig_lung_nasal_day5_ds = cs.create_pig_lung_nasal_day5()


# Fig1.fig1_develop(guinea_pig_nasal_ds, ferret_nasal_ds, pig_nasal_ds)
# Fig2f.fig2_develop(ferret_lung_nasal_ds)
# Fig3.fig3_develop(ferret_lung_nasal_ds)
# Fig4f.fig4_develop(pig_lung_day3_ds, pig_lung_nasal_day5_ds)
# Fig5.fig5_develop(pig_lung_day3_ds, pig_lung_nasal_day5_ds)

# Fig1sup.fig1sup_develop(guinea_pig_nasal_ds, ferret_nasal_ds, pig_nasal_ds)
# Fig2sup.fig2sup_develop(pig_lung_day3_ds, pig_lung_nasal_day5_ds, ferret_lung_nasal_ds)
# Fig3sup.fig3sup_develop()
