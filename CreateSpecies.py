import SpeciesData as sd
"""Creates species as pickle is not working"""


def create_ferret_nasal(is_connectivity_index=False):
    file = 'Ferret Nasal Wash & Lung NL09 WT_VAR Genotype Tables.xls'
    sheet_number = 0
    species_name = 'Ferret Nasal'
    x_dimension_starts = [[1, 10, 19]] * 10
    x_dimension_finishes = [[8, 17, 26]] * 10
    y_dimension_starts = [[3] * 3, [26] * 3, [49] * 3, [72] * 3, [95] * 3, [118] * 3, [143] * 3, [166] * 3, [189] * 3,
                          [212] * 3]
    y_dimension_finishes = [[23] * 3, [46] * 3, [69] * 3, [92] * 3, [115] * 3, [138] * 3, [163] * 3, [186] * 3,
                            [209] * 3,
                            [232] * 3]
    individual_names = [f'F{number + 1}' for number in range(10)]
    data_points = [['Day 1', 'Day 3', 'Day 5']] * 10
    ferret_nasal_ds = sd.Species(file, sheet_number, species_name, x_dimension_starts, x_dimension_finishes,
                                 y_dimension_starts, y_dimension_finishes, individual_names, data_points, split_point=6,
                                 is_connectivity_index=is_connectivity_index)
    return ferret_nasal_ds


def create_ferret_lung(is_connectivity_index=False):
    file = 'Ferret Nasal Wash & Lung NL09 WT_VAR Genotype Tables.xls'
    sheet_number = 1
    species_name = 'Ferret Lung'
    x_dimension_starts = [[1], [1], [10], [10], [19], [19], [28], [28], [1], [1], [10], [10], [19], [28]]
    x_dimension_finishes = [[8], [8], [17], [17], [26], [26], [35], [35], [8], [8], [17], [17], [26], [35]]
    y_dimension_starts = [[2], [25], [2], [25], [2], [25], [2], [25], [48], [71], [48], [71], [48], [48]]
    y_dimension_finishes = [[22], [45], [22], [45], [22], [45], [22], [45], [68], [91], [68], [91], [68], [68]]
    individual_names = ['F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'F17', 'F18', 'F19', 'F20', 'F21', 'F22', 'F23',
                        'F24']
    # individual_names = ['F12', 'F13', 'F14', 'F16', 'F17', 'F19', 'F20', 'F22', 'F23', 'F25', 'F27', 'F28', 'F30',
    #                     'F34']
    data_points = [['Day 1']] * 2 + [['Day 2']] * 2 + [['Day 3']] * 2 + [['Day 4']] * 2 + [['Day 1']] * 2 + \
                  [['Day 2']] * 2 + [['Day 3']] + [['Day 4']]
    ferret_lung_ds = sd.Species(file, sheet_number, species_name, x_dimension_starts, x_dimension_finishes,
                                y_dimension_starts, y_dimension_finishes, individual_names, data_points, split_point=8,
                                is_connectivity_index=is_connectivity_index)
    return ferret_lung_ds


def create_ferret_lung_nasal(is_connectivity_index=False):
    file = 'Ferret Nasal Turbinate & Lung NL09 WT_VAR Genotype Tables.xls'
    sheet_number = 0
    species_name = 'Ferret Lung and Nasal Turbinate'
    x_dimension_starts = [[1, 1], [10, 10], [19, 19], [28, 28], [37, 37], [46, 46], [55, 55], [64, 64], [73, 73], [82, 82], [91, 91], [100, 100], [109, 109], [118, 118]]
    x_dimension_finishes = [[8, 8], [17, 17], [26, 26], [35, 35], [44, 44], [53, 53], [62, 62], [71, 71], [80, 80], [89, 89], [98, 98], [107, 107], [116, 116], [125, 125]]
    y_dimension_starts = [[2, 25]]*14
    y_dimension_finishes = [[22, 45]]*14
    individual_names = ['F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'F17', 'F18', 'F19', 'F20', 'F21', 'F22', 'F23',
                        'F24']
    # individual_names = ['F12', 'F13', 'F14', 'F16', 'F17', 'F19', 'F20', 'F22', 'F23', 'F25', 'F27', 'F28', 'F30',
    #                     'F34']
    data_points = [['Lung', 'Nasal Turbinate']]*14
    ferret_lung_nasal_ds = sd.Species(file, sheet_number, species_name, x_dimension_starts, x_dimension_finishes,
                                      y_dimension_starts, y_dimension_finishes, individual_names, data_points,
                                      split_point=8, is_connectivity_index=is_connectivity_index)
    return ferret_lung_nasal_ds


def create_guinea_pig_nasal(is_connectivity_index=False):
    file = 'Guinea Pig Nasal Wash Reassortment Genotype Tables.xls'
    sheet_number = 0
    species_name = 'Guinea Pig Nasal'
    x_dimension_starts = [[1, 10, 19]] * 10
    x_dimension_finishes = [[8, 17, 26]] * 10
    y_dimension_starts = [[3] * 3, [26] * 3, [49] * 3, [72] * 3, [95] * 3, [118] * 3, [145] * 3, [167] * 3, [191] * 3,
                          [214] * 3]
    y_dimension_finishes = [[23, 22, 23], [45, 46, 46], [69] * 3, [91, 92, 91], [115, 114, 113], [138, 137, 138],
                            [162, 162, 164], [187] * 3, [210, 211, 211], [233, 231, 234]]
    individual_names = [f'GP{number + 1}' for number in range(10)]
    data_points = [['Day 1', 'Day 2', 'Day 4']] * 6 + [['Day 1', 'Day 2', 'Day 3']] * 4
    guinea_pig_nasal_ds = sd.Species(file, sheet_number, species_name, x_dimension_starts, x_dimension_finishes,
                                     y_dimension_starts, y_dimension_finishes, individual_names, data_points,
                                     split_point=6, is_connectivity_index=is_connectivity_index)
    return guinea_pig_nasal_ds


def create_pig_lung_day3(is_connectivity_index=False):
    file = 'Pig Lung Homogenate Genotype Tables.xls'
    sheet_number = 0
    species_name = 'Pig Lung Day 3'
    x_dimension_starts = [[1] * 7, [10] * 7, [19] * 7]
    x_dimension_finishes = [[8] * 7, [17] * 7, [26] * 7]
    y_dimension_starts = [[2, 25, 48, 71, 94, 117, 140]] * 3
    y_dimension_finishes = [[22, 45, 68, 91, 114, 137, 160]] * 3
    individual_names = ['P1', 'P2', 'P3']
    # individual_names = ['P79', 'P84', 'P73']
    data_points = [['L.Apical', 'R.Apical', 'L.Cardiac', 'R.Cardiac', 'L.Diaphr', 'R.Diaphr', 'Intermed']] * 3
    pig_lung_day3_ds = sd.Species(file, sheet_number, species_name, x_dimension_starts, x_dimension_finishes,
                                  y_dimension_starts, y_dimension_finishes, individual_names, data_points,
                                  is_spatial=True, is_connectivity_index=is_connectivity_index)
    return pig_lung_day3_ds


def create_pig_lung_day5(is_connectivity_index=False):
    file = 'Pig Lung Homogenate Genotype Tables.xls'
    sheet_number = 1
    species_name = 'Pig Lung Day 5'
    x_dimension_starts = [[1] * 7, [10] * 7, [19] * 7]
    x_dimension_finishes = [[8] * 7, [17] * 7, [26] * 7]
    y_dimension_starts = [[2, 25, 94, 117, 48, 71, 140]] * 3
    y_dimension_finishes = [[22, 45, 114, 137, 68, 91, 160]] * 3
    individual_names = ['P4', 'P5', 'P6']
    # individual_names = ['P55', 'P52', 'P76']
    # data_points = [['L.Apical', 'R.Apical', 'L.Diaphr', 'R.Diaphr', 'L.Cardiac', 'R.Cardiac', 'Intermed']] * 3
    data_points = [['L.Apical', 'R.Apical', 'L.Cardiac', 'R.Cardiac', 'L.Diaphr', 'R.Diaphr', 'Intermed']] * 3
    pig_lung_day5_ds = sd.Species(file, sheet_number, species_name, x_dimension_starts, x_dimension_finishes,
                                  y_dimension_starts, y_dimension_finishes, individual_names, data_points,
                                  is_spatial=True, is_connectivity_index=is_connectivity_index)
    return pig_lung_day5_ds


def create_pig_nasal(is_connectivity_index=False):
    file = 'Pig Nasal Swab Genotype Tables.xls'
    sheet_number = 0
    species_name = 'Pig Nasal'
    x_dimension_starts = [[1, 10, 19]] * 6
    x_dimension_finishes = [[8, 17, 26]] * 6
    y_dimension_starts = [[25] * 3, [2] * 3, [71] * 3, [48] * 3, [94] * 3, [117] * 3]
    y_dimension_finishes = [[45] * 3, [22] * 3, [91] * 3, [68] * 3, [114] * 3, [137, 137, 135]]
    individual_names = ['P4', 'P5', 'P6', 'P7', 'P8', 'P9']
    # individual_names = ['P55', 'P52', 'P76', 'P66', 'P88', 'P89']
    data_points = [['Day 1', 'Day 3', 'Day 5']] * 6
    pig_nasal_ds = sd.Species(file, sheet_number, species_name, x_dimension_starts, x_dimension_finishes,
                              y_dimension_starts, y_dimension_finishes, individual_names, data_points,
                              is_connectivity_index=is_connectivity_index)
    return pig_nasal_ds


def create_pig_lung_nasal_day5(is_connectivity_index=False):
    file = 'Pig Lung Homogenate Genotype Tables.xls'
    sheet_number = 1
    species_name = 'Pig Lung and Nasal Day 5'
    x_dimension_starts = [[1] * 8, [10] * 8, [19] * 8]
    x_dimension_finishes = [[8] * 8, [17] * 8, [26] * 8]
    y_dimension_starts = [[2, 25, 94, 117, 48, 71, 140, 163]] * 3
    y_dimension_finishes = [[22, 45, 114, 137, 68, 91, 160, 183]] * 3
    individual_names = ['P4', 'P5', 'P6']
    # individual_names = ['P55', 'P52', 'P76']
    # data_points = [['L.Apical', 'R.Apical', 'L.Diaphr', 'R.Diaphr', 'L.Cardiac', 'R.Cardiac', 'Intermed', 'Nasal']] * 3
    data_points = [['L.Apical', 'R.Apical', 'L.Cardiac', 'R.Cardiac', 'L.Diaphr', 'R.Diaphr', 'Intermed', 'Nasal']] * 3
    pig_lung_day5_ds = sd.Species(file, sheet_number, species_name, x_dimension_starts, x_dimension_finishes,
                                  y_dimension_starts, y_dimension_finishes, individual_names, data_points,
                                  is_spatial=True, is_connectivity_index=is_connectivity_index)
    return pig_lung_day5_ds


def create_quail_tracheal(is_connectivity_index=False):
    file = 'Quail Tracheal Swab Genotype Tables.xls'
    sheet_number = 0
    species_name = 'Quail Tracheal'
    x_dimension_starts = [[1, 10, 19]] * 3 + [[1, 10]] + [[1, 10, 19]] * 6
    x_dimension_finishes = [[8, 17, 26]] * 3 + [[8, 17]] + [[8, 17, 26]] * 6
    y_dimension_starts = [[3] * 3, [26] * 3, [49] * 3, [72] * 2, [95] * 3, [118] * 3, [141] * 3, [164] * 3, [187] * 3,
                          [209] * 3]
    y_dimension_finishes = [[23, 22, 23], [45, 45, 46], [69, 69, 68], [92] * 2, [115, 112, 114], [138] * 3, [161] * 3,
                            [184] * 3, [205, 205, 206], [229, 230, 230]]
    individual_names = ['Q19', 'Q20', 'Q21', 'Q22', 'Q23', 'Q24', 'Q41', 'Q42', 'Q44', 'Q45']
    data_points = [['Day 1', 'Day 3', 'Day 5']] * 3 + [['Day 1', 'Day 3']] + [['Day 1', 'Day 3', 'Day 5']] * 6
    quail_tracheal_ds = sd.Species(file, sheet_number, species_name, x_dimension_starts, x_dimension_finishes,
                                   y_dimension_starts, y_dimension_finishes, individual_names, data_points,
                                   split_point=6, is_connectivity_index=is_connectivity_index)
    return quail_tracheal_ds


def create_guinea_pig_verify():
    """Creates guinea pig species from previous study to verify code is working"""
    file = 'GP Reassortment_H3N8_H4N6.xls'
    sheet_number = 0
    species_name = 'Guinea Pig'
    x_dimension_starts = [[1, 10, 19]] * 8
    x_dimension_finishes = [[8, 17, 26]] * 8
    y_dimension_starts = [[2] * 3, [25] * 3, [48] * 3, [71] * 3, [94] * 3, [117] * 3, [140] * 3, [163] * 3]
    y_dimension_finishes = [[22] * 3, [45] * 3, [68] * 3, [91] * 3, [114] * 3, [137] * 3, [160] * 3, [183] * 3]
    individual_names = [f'GP#{number + 1}' for number in range(8)]
    data_points = [['Day 1', 'Day 2', 'Day 3']] * 8
    GP_ds = sd.Species(file, sheet_number, species_name, x_dimension_starts, x_dimension_finishes, y_dimension_starts,
                       y_dimension_finishes, individual_names, data_points)
    return GP_ds
