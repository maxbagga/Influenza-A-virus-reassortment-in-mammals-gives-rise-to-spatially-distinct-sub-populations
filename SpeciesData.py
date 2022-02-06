import xlrd
import IndividualData as id
import statistics
import pandas as pd
"""SpeciesData Creates Species"""


class Species:
    def __init__(self, file, sheet_number, species_name, x_dimension_starts, x_dimension_finishes, y_dimension_starts,
                 y_dimension_finishes, individual_names, data_points, is_spatial=False, split_point=None,
                 is_connectivity_index=False, is_combined_sheets=False, ds_split_point=None):
        """Species defines an entire sheet Data"""
        self.file = file
        self.sheet_number = sheet_number
        self.wb = xlrd.open_workbook(self.file)
        if not is_combined_sheets:
            self.sheet = self.wb.sheet_by_index(self.sheet_number)
        else:
            sheets = []
            for sheet in sheet_number:
                sheets.append(self.wb.sheet_by_index(sheet))
            self.sheet = sheets
        self.species_name = species_name
        self.x_dimension_starts = x_dimension_starts
        self.x_dimension_finishes = x_dimension_finishes
        self.y_dimension_starts = y_dimension_starts
        self.y_dimension_finishes = y_dimension_finishes
        self.individual_names = individual_names
        self.number_of_individuals = len(individual_names)
        self.data_points = data_points
        self.is_connectivity_index = is_connectivity_index
        self.is_combined_sheets = is_combined_sheets
        self.ds_split_point = ds_split_point
        self.individuals = self.get_individuals()
        self.len_data_points = [len(data_p) for data_p in self.data_points]
        self.max_len_data_points = max(self.len_data_points)
        self.max_len_data_points_location = self.get_max_len_location()
        self.is_spatial = is_spatial
        self.split_point = split_point
        # self.plaques_string, self.plaques_float, self.unique_genotypes_df = self.get_plaques()
        # self.number_of_plaques = len(self.plaques_string)
        if self.is_connectivity_index:
            self.connectivity_indices, self.connectivity_index_graphs, self.max_connectivity_index,\
                self.max_connectivity_index_graph = self.get_connectivity_index()
            self.average_connectivity_index = statistics.mean(self.connectivity_indices)
        self.ten_squared, self.ten_fifth = self.get_split()

    def get_individuals(self):
        """Passes info and creates individuals in species"""
        individual_list = []
        for individual in range(self.number_of_individuals):
            if not self.is_combined_sheets:
                sheet = self.sheet
            else:
                if individual < self.ds_split_point:
                    sheet = self.sheet[0]
                else:
                    sheet = self.sheet[1]
            individual_list.append(id.Individual(self.wb, sheet, self.individual_names[individual],
                                                 self.x_dimension_starts[individual],
                                                 self.x_dimension_finishes[individual],
                                                 self.y_dimension_starts[individual],
                                                 self.y_dimension_finishes[individual], self.data_points[individual],
                                                 is_connectivity_index=self.is_connectivity_index,
                                                 is_combined_sheets=self.is_combined_sheets))
        return individual_list

    def get_max_len_location(self):
        """Gets first location of individual with the maximum number of datapoints in species"""
        for index in range(self.number_of_individuals):
            if len(self.data_points[index]) == self.max_len_data_points:
                return index

    def get_split(self):
        """Splits the data into different inoculation doses"""
        if self.split_point is not None:
            ten_squared = Species(self.file, self.sheet_number, self.species_name,
                                  self.x_dimension_starts[:self.split_point],
                                  self.x_dimension_finishes[:self.split_point],
                                  self.y_dimension_starts[:self.split_point],
                                  self.y_dimension_finishes[:self.split_point],
                                  self.individual_names[:self.split_point], self.data_points[:self.split_point],
                                  is_spatial=self.is_spatial)
            ten_fifth = Species(self.file, self.sheet_number, self.species_name,
                                self.x_dimension_starts[self.split_point:],
                                self.x_dimension_finishes[self.split_point:],
                                self.y_dimension_starts[self.split_point:],
                                self.y_dimension_finishes[self.split_point:],
                                self.individual_names[self.split_point:], self.data_points[self.split_point:],
                                is_spatial=self.is_spatial)
            return ten_squared, ten_fifth
        else:
            return None, None

    def get_connectivity_index(self):
        """Compiles the optimal connectivity indices & graph for individuals and the maximum values for both"""
        connectivity_indices = []
        connectivity_index_graphs = []
        max_connectivity_index = 0
        max_connectivity_index_graph = None
        for individual in self.individuals:
            connectivity_indices.append(individual.optimal_connectivity_index)
            connectivity_index_graphs.append(individual.optimal_connectivity_index_graph)
            if individual.optimal_connectivity_index > max_connectivity_index:
                max_connectivity_index = individual.optimal_connectivity_index
                max_connectivity_index_graph = individual.optimal_connectivity_index_graph.__class__()
                max_connectivity_index_graph.add_nodes_from(individual.optimal_connectivity_index_graph)
                max_connectivity_index_graph.add_edges_from(individual.optimal_connectivity_index_graph.edges)
        return connectivity_indices, connectivity_index_graphs, max_connectivity_index, max_connectivity_index_graph

    def get_plaques(self):
        """Combines plaques from each individual and also forms data frame with ids and appearances"""
        plaques_string = []
        plaques_float = []
        gas = [individual.genotype_analyzer for individual in self.individuals]
        for ga1 in range(self.number_of_individuals):
            for ga2 in range(self.number_of_individuals):
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
        for individual in self.individuals:
            plaques_string += individual.plaques_string
            plaques_float += individual.plaques_float
        return plaques_string, plaques_float, unique_genotypes_df

    def __add__(self, other):
        new_ds = Species(self.file, [self.sheet_number, other.sheet_number], [self.species_name, other.species_name],
                         self.x_dimension_starts + other.x_dimension_starts, self.x_dimension_finishes +
                         other.x_dimension_finishes, self.y_dimension_starts + other.y_dimension_starts,
                         self.y_dimension_finishes + other.y_dimension_finishes, self.individual_names +
                         other.individual_names, self.data_points + other.data_points, is_spatial=self.is_spatial,
                         split_point=self.split_point, is_connectivity_index=self.is_connectivity_index,
                         is_combined_sheets=True, ds_split_point=self.number_of_individuals)
        return new_ds
