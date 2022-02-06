import matplotlib.pyplot as plt
import Fig2
import Fig6
from scipy.stats import ttest_rel
import pandas as pd


def fig1_develop(gp_ds, f_ds, p_ds):
    fig, axs = plt.subplots(nrows=4, ncols=4, figsize=(40, 40))
    writer = pd.ExcelWriter(f'/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 '
                            f'Results/Figure Data/Fig 1 Data.xlsx', engine='openpyxl')
    gp_pf_all, gp_pf_1, gp_pf_2, gp_richness_all, gp_richness_1, gp_richness_2, gp_diversity_all, gp_diversity_1, \
    gp_diversity_2, = Fig2.fig2_develop(gp_ds, axs=axs[0, 0:4], target=9, labels=['A', 'B', 'C', 'D'], code=0,
                                        tick_points=5, plot_points=[1, 2, 4], split_lines=True, excel_writer=writer)
    f_pf, f_richness, f_diversity = Fig2.fig2_develop(f_ds, axs=axs[1, 0:4], target=0, labels=['E', 'F', 'G', 'H'],
                                                      code=1, tick_points=5, plot_points=[1, 3, 5], excel_writer=writer)
    p_pf, p_richness, p_diversity = Fig2.fig2_develop(p_ds, axs=axs[2, 0:4], target=4, labels=['I', 'J', 'K', 'L'],
                                                      code=2, tick_points=5, plot_points=[1, 3, 5], is_x_label=True,
                                                      excel_writer=writer)
    #Figure data
    gp_pf_all.to_excel(writer, sheet_name='Fig 1B')
    gp_richness_all.to_excel(writer, sheet_name='Fig 1C')
    gp_diversity_all.to_excel(writer, sheet_name='Fig 1D')
    f_pf.to_excel(writer, sheet_name='Fig 1F')
    f_richness.to_excel(writer, sheet_name='Fig 1G')
    f_diversity.to_excel(writer, sheet_name='Fig 1H')
    p_pf.to_excel(writer, sheet_name='Fig 1J')
    p_richness.to_excel(writer, sheet_name='Fig 1K')
    p_diversity.to_excel(writer, sheet_name='Fig 1L')

    gp_pf_all = gp_pf_all.values.tolist()
    gp_pf_1 = gp_pf_1.values.tolist()
    gp_pf_2 = gp_pf_2.values.tolist()

    gp_richness_all = gp_richness_all.values.tolist()
    gp_richness_1 = gp_richness_1.values.tolist()
    gp_richness_2 = gp_richness_2.values.tolist()

    gp_diversity_all = gp_diversity_all.values.tolist()
    gp_diversity_1 = gp_diversity_1.values.tolist()
    gp_diversity_2 = gp_diversity_2.values.tolist()

    f_pf = f_pf.values.tolist()
    f_richness = f_richness.values.tolist()
    f_diversity = f_diversity.values.tolist()

    p_pf = p_pf.values.tolist()
    p_richness = p_richness.values.tolist()
    p_diversity = p_diversity.values.tolist()

    print('Guinea Pig paired t-tests-')
    print(f'Parental Frequency Day 1 Day 2: {ttest_rel(gp_pf_all[0], gp_pf_all[1])[1]}')
    print(f'Parental Frequency Day 1 Day 3: {ttest_rel(gp_pf_2[0], gp_pf_2[2])[1]}')
    print(f'Parental Frequency Day 1 Day 4: {ttest_rel(gp_pf_1[0], gp_pf_1[2])[1]}')
    print(f'Parental Frequency Day 2 Day 3: {ttest_rel(gp_pf_2[1], gp_pf_2[2])[1]}')
    print(f'Parental Frequency Day 2 Day 4: {ttest_rel(gp_pf_1[1], gp_pf_1[2])[1]}')

    print(f'Richness Day 1 Day 2: {ttest_rel(gp_richness_all[0], gp_richness_all[1])[1]}')
    print(f'Richness Day 1 Day 3: {ttest_rel(gp_richness_2[0], gp_richness_2[2])[1]}')
    print(f'Richness Day 1 Day 4: {ttest_rel(gp_richness_1[0], gp_richness_1[2])[1]}')
    print(f'Richness Day 2 Day 3: {ttest_rel(gp_richness_2[1], gp_richness_2[2])[1]}')
    print(f'Richness Day 2 Day 4: {ttest_rel(gp_richness_1[1], gp_richness_1[2])[1]}')

    print(f'Diversity Day 1 Day 2: {ttest_rel(gp_diversity_all[0], gp_diversity_all[1])[1]}')
    print(f'Diversity Day 1 Day 3: {ttest_rel(gp_diversity_2[0], gp_diversity_2[2])[1]}')
    print(f'Diversity Day 1 Day 4: {ttest_rel(gp_diversity_1[0], gp_diversity_1[2])[1]}')
    print(f'Diversity Day 2 Day 3: {ttest_rel(gp_diversity_2[1], gp_diversity_2[2])[1]}')
    print(f'Diversity Day 2 Day 4: {ttest_rel(gp_diversity_1[1], gp_diversity_1[2])[1]}')

    print('Ferret paired t-tests-')
    print(f'Parental Frequency Day 1 Day 3: {ttest_rel(f_pf[0], f_pf[1])[1]}')
    print(f'Parental Frequency Day 1 Day 5: {ttest_rel(f_pf[0], f_pf[2])[1]}')
    print(f'Parental Frequency Day 3 Day 5: {ttest_rel(f_pf[1], f_pf[2])[1]}')

    print(f'Richness Day 1 Day 3: {ttest_rel(f_richness[0], f_richness[1])[1]}')
    print(f'Richness Day 1 Day 5: {ttest_rel(f_richness[0], f_richness[2])[1]}')
    print(f'Richness Day 3 Day 5: {ttest_rel(f_richness[1], f_richness[2])[1]}')

    print(f'Diversity Day 1 Day 3: {ttest_rel(f_diversity[0], f_diversity[1])[1]}')
    print(f'Diversity Day 1 Day 5: {ttest_rel(f_diversity[0], f_diversity[2])[1]}')
    print(f'Diversity Day 3 Day 5: {ttest_rel(f_diversity[1], f_diversity[2])[1]}')

    print('Pig paired t-tests-')
    print(f'Parental Frequency Day 1 Day 3: {ttest_rel(p_pf[0], p_pf[1])[1]}')
    print(f'Parental Frequency Day 1 Day 5: {ttest_rel(p_pf[0], p_pf[2])[1]}')
    print(f'Parental Frequency Day 3 Day 5: {ttest_rel(p_pf[1], p_pf[2])[1]}')

    print(f'Richness Day 1 Day 3: {ttest_rel(p_richness[0], p_richness[1])[1]}')
    print(f'Richness Day 1 Day 5: {ttest_rel(p_richness[0], p_richness[2])[1]}')
    print(f'Richness Day 3 Day 5: {ttest_rel(p_richness[1], p_richness[2])[1]}')

    print(f'Diversity Day 1 Day 3: {ttest_rel(p_diversity[0], p_diversity[1])[1]}')
    print(f'Diversity Day 1 Day 5: {ttest_rel(p_diversity[0], p_diversity[2])[1]}')
    print(f'Diversity Day 3 Day 5: {ttest_rel(p_diversity[1], p_diversity[2])[1]}')

    Fig6.fig6_develop(gp_ds, f_ds, p_ds, axs=(axs[3, 1], axs[3, 2], axs[3, 3]), excel_writer=writer)
    writer.save()

    axs[3, 0].axis('off')

    plt.savefig('/Users/maxbagga/OneDrive - Emory University/Box/Mallard Phylogeny Data/Manuscript #2 Results/'
                'Fig 1.tiff', dpi=300, format='tiff')
