import IndividualDayData as idd
from itertools import combinations
import pandas as pd
import networkx as NX
import LDObject as ldo
import UniqueGenotype as ug
import statistics
import numpy as np
import matplotlib.pyplot as plt
import random
"""IndividualData creates Individual"""


class Individual:
    def __init__(self, wb, sheet, name, x_dimension_starts, x_dimension_finishes, y_dimension_starts,
                 y_dimension_finishes, data_points, is_connectivity_index=False, is_combined_sheets=False,
                 split_sheets=False, is_simulation=False, override_plaques=None):
        """Individual defines one row in genotype table"""
        self.wb = wb
        self.sheet = sheet
        self.name = name
        self.x_dimension_starts = x_dimension_starts
        self.x_dimension_finishes = x_dimension_finishes
        self.y_dimension_starts = y_dimension_starts
        self.y_dimension_finishes = y_dimension_finishes
        self.data_points = data_points
        self.number_of_days = len(data_points)
        self.is_combined_sheets = is_combined_sheets
        self.spit_sheets = split_sheets
        self.is_simulation = is_simulation
        self.override_plaques = override_plaques
        self.days = self.get_days(override_plaques=override_plaques)
        self.plaques_string, self.plaques_float, self.unique_genotypes_df = self.get_plaques()
        self.number_of_plaques = len(self.plaques_string)
        self.ldos = self.get_ldos()
        self.unique_genotypes, self.appearances = self.get_unique_genotypes_appearances()
        self.number_of_unique_genotypes = len(self.unique_genotypes)
        self.unique_genotype_objects = [ug.Genotype(self.unique_genotypes[genotype], self.appearances[genotype],
                                                    self.ldos, self.number_of_plaques) for genotype in
                                        range(self.number_of_unique_genotypes)]
        self.genotype_analyzer = ug.GenotypeAnalyzer(self.unique_genotype_objects)
        self.q_0 = self.genotype_analyzer.unique
        self.q_1 = self.genotype_analyzer.shannon_wiener
        self.q_infinity = 1/(max(self.appearances)/self.number_of_plaques)
        self.beta_diversity_0, self.beta_diversity_1, self.beta_diversity_infinity, self.Meff_0_standardized, \
            self.Meff_1_standardized, self.Meff_infinity_standardized = self.get_beta_diversities()
        self.is_connectivity_index = is_connectivity_index
        if self.is_connectivity_index:
            self.optimal_connectivity_index_graph, self.optimal_connectivity_index = \
                self.get_optimal_connectivity_index_graph()
        self.ug_intersection = self.get_ug_intersection()
        self.jaccard_index = len(self.ug_intersection)/self.number_of_unique_genotypes

    def get_days(self, override_plaques=None):
        """Passes info and creates individual days in individual"""
        day_list = []
        for day in range(self.number_of_days):
            if not self.spit_sheets:
                sheet = self.sheet
            else:
                sheet = self.sheet[day]
            if override_plaques is None:
                day_list.append(idd.IndividualDay(self.wb, sheet, self.name, self.x_dimension_starts[day],
                                                  self.x_dimension_finishes[day], self.y_dimension_starts[day],
                                                  self.y_dimension_finishes[day], self.data_points[day],
                                                  is_simulation=self.is_simulation))
            else:
                day_list.append(idd.IndividualDay(self.wb, sheet, self.name, self.x_dimension_starts[day],
                                                  self.x_dimension_finishes[day], self.y_dimension_starts[day],
                                                  self.y_dimension_finishes[day], self.data_points[day],
                                                  is_simulation=self.is_simulation,
                                                  override_plaques=override_plaques[day]))
        return day_list

    def get_plaques(self):
        """Combines plaques from each day in individual and also forms data frame with ids and appearances"""
        plaques_string = []
        plaques_float = []
        gas = [day.genotype_analyzer for day in self.days]
        for ga1 in range(self.number_of_days):
            for ga2 in range(self.number_of_days):
                if ga1 != ga2:
                    gas[ga1], gas[ga2] = gas[ga1] + gas[ga2]
        ids = gas[0].ids
        appearances = []
        for id in range(len(ids)):
            appearance = 0
            for ga in gas:
                appearance += ga.appearances[id]
            appearances.append(appearance)
        unique_genotypes_df = pd.DataFrame({'Unique Genotypes': ids,
                                            'Appearances': appearances})
        for day in self.days:
            plaques_string += day.plaques_string
            plaques_float += day.plaques_float
        return plaques_string, plaques_float, unique_genotypes_df

    def get_optimal_connectivity_index_graph(self):
        """Return optimal graph and connectivity for individual based on having 4 similarities out of 6"""
        optimal_graph = None
        optimal_score = 0
        for combination1 in list(combinations(range(8), 6)):
            g = NX.Graph()
            g.add_nodes_from(self.unique_genotypes_df['Unique Genotypes'])
            for combination2 in list(combinations(self.unique_genotypes_df['Unique Genotypes'], 2)):
                count = 0
                for index in combination1:
                    if combination2[0][index] == combination2[1][index]:
                        count += 1
                if count >= 4:
                    g.add_edge(combination2[0], combination2[1])
                    g.edges[combination2[0], combination2[1]]['Appearances'] = \
                        self.unique_genotypes_df[self.unique_genotypes_df['Unique Genotypes'] ==
                                                 combination2[0]]['Appearances'].values[0] + \
                        self.unique_genotypes_df[self.unique_genotypes_df['Unique Genotypes'] ==
                                                 combination2[1]]['Appearances'].values[0]
            connectivity_index = 0
            for edge in g.edges(data=True):
                degrees = g.degree
                degree1 = 0
                degree2 = 0
                for degree in degrees:
                    if edge[0] == degree[0]:
                        degree1 = degree[1]
                    elif edge[1] == degree[0]:
                        degree2 = degree[1]
                weight = edge[2].get('Appearances')
                connectivity_index += weight/((degree1 * degree2) ** .5)
            if connectivity_index > optimal_score:
                optimal_score = connectivity_index
                optimal_graph = g.__class__()
                optimal_graph.add_nodes_from(g)
                optimal_graph.add_edges_from(g.edges)
        return optimal_graph, optimal_score

    def get_unique_genotypes_appearances(self):
        """Identifies the unique genotypes and their appearances and returns them"""
        unique_genotypes = []
        for plaque in self.plaques_string:
            if plaque in unique_genotypes:
                pass
            else:
                unique_genotypes.append(plaque)
        appearances = []
        for plaque in unique_genotypes:
            appearances.append(self.plaques_string.count(plaque))
        return unique_genotypes, appearances

    def get_ldos(self):
        """Sorts plaques_float by segment which is used in certain statistics"""
        segment_values_list = []
        for segment in range(8):
            current_segment = []
            for plaque in range(self.number_of_plaques):
                current_segment.append(self.plaques_float[plaque][segment])
            segment_values_list.append(current_segment)
        ldos = []
        for segment in segment_values_list:
            ldos.append(ldo.LDObject(segment))
        return ldos

    def get_beta_diversities(self):
        """Returns Beta Diversity Values"""
        beta_diversity_0 = self.q_0/statistics.mean([day.q_0 for day in self.days]) if \
            self.q_0/statistics.mean([day.q_0 for day in self.days]) >= 1 else 1
        beta_diversity_1 = self.q_1/statistics.mean([day.q_1 for day in self.days]) if \
            self.q_1/statistics.mean([day.q_1 for day in self.days]) >= 1 else 1
        beta_diversity_infinity = self.q_infinity/statistics.mean([day.q_infinity for day in self.days]) if \
            self.q_infinity/statistics.mean([day.q_infinity for day in self.days]) >= 1 else 1
        Meff_0_standardized = (beta_diversity_0 - 1) / (self.number_of_days - 1) if \
            (beta_diversity_0 - 1) / (self.number_of_days - 1) <= 1 else 1
        Meff_1_standardized = (beta_diversity_1 - 1) / (self.number_of_days - 1) if \
            (beta_diversity_1 - 1) / (self.number_of_days - 1) <= 1 else 1
        Meff_infinity_standardized = (beta_diversity_infinity - 1) / (self.number_of_days - 1) if \
            (beta_diversity_infinity - 1) / (self.number_of_days - 1) <= 1 else 1
        return beta_diversity_0, beta_diversity_1, beta_diversity_infinity, Meff_0_standardized, Meff_1_standardized, \
               Meff_infinity_standardized

    def get_ug_intersection(self):
        """Returns number of unique genotypes in union and intersection of days"""
        intersection_of_ugs = [ug for ug in self.unique_genotypes]
        for day in self.days:
            temp_intersection_of_ugs = [ug for ug in intersection_of_ugs]
            for ug1 in intersection_of_ugs:
                if ug1 not in day.unique_genotypes:
                    temp_intersection_of_ugs.remove(ug1)
            intersection_of_ugs = [ug for ug in temp_intersection_of_ugs]
        return intersection_of_ugs

    def get_clustering_matrix(self):
        uh_matrix = np.zeros((2, 2, 2, 2, 2, 2, 2, 2))
        for plaque in self.plaques_float:
            current_plaque = [int(segment) for segment in plaque]
            uh_matrix[current_plaque[0]][current_plaque[1]][current_plaque[2]][current_plaque[3]][current_plaque[4]][current_plaque[5]][current_plaque[6]][current_plaque[7]] += 1

    def reshuffle_plaques(self):
        """Shuffles plaques across genotype tables in order to test compartmentalization"""
        plaques_string_reshuffled = random.sample(self.plaques_string, len(self.plaques_string))
        count = 0
        day_plaques_string_reshuffled = []
        for day in range(self.number_of_days):
            end_index = self.y_dimension_finishes[day] - self.y_dimension_starts[day] + 1
            day_plaques_string_reshuffled.append(plaques_string_reshuffled[count:count + end_index])
            count += end_index
        return self.get_days(override_plaques=day_plaques_string_reshuffled)






