import LDObject as ldo
import UniqueGenotype as ug
import pandas as pd
import random
import numpy as np

"""IndividualDayData starts data analysis process"""


plaque_generation_size = 1000


class IndividualDay:
    def __init__(self, wb, sheet, name, x_dimension_start, x_dimension_finish, y_dimension_start, y_dimension_finish,
                 day, is_simulation=False, override_plaques=None):
        """Individual defines one table"""
        self.wb = wb
        self.sheet = sheet
        self.name = name
        self.x_dimension_start = x_dimension_start
        self.x_dimension_finish = x_dimension_finish
        self.y_dimension_start = y_dimension_start
        self.y_dimension_finish = y_dimension_finish
        self.day = day
        self.is_simulation = is_simulation
        if override_plaques is None:
            self.plaques_string, self.plaques_float = self.get_plaques()
        else:
            self.plaques_string = override_plaques
            self.plaques_float = [[float(letter) for letter in plaque] for plaque in self.plaques_string]
        self.number_of_plaques = len(self.plaques_string)
        self.ldos = self.get_ldos()
        self.full_name = f'{self.name} {self.day}'
        self.unique_genotypes, self.appearances = self.get_unique_genotypes_appearances()
        self.unique_genotypes_df = pd.DataFrame({'Unique Genotypes': self.unique_genotypes,
                                                 'Appearances': self.appearances})
        self.number_of_unique_genotypes = len(self.unique_genotypes)
        self.unique_genotype_objects = [ug.Genotype(self.unique_genotypes[genotype], self.appearances[genotype],
                                                    self.ldos, self.number_of_plaques) for genotype in
                                        range(self.number_of_unique_genotypes)]
        self.genotype_analyzer = ug.GenotypeAnalyzer(self.unique_genotype_objects)
        self.q_0 = self.genotype_analyzer.unique
        self.q_1 = self.genotype_analyzer.shannon_wiener
        self.q_infinity = 1/(max(self.appearances)/self.number_of_plaques)

    def get_plaques(self):
        """Process spreadsheet and creates both string and float plaques"""
        plaque_string_list = []
        plaque_float_list = []
        for y in range(self.y_dimension_finish - self.y_dimension_start+1):
            cell_string = ''
            cell_float = []
            go_on = True
            for x in range(self.x_dimension_finish - self.x_dimension_start + 1):
                # print(y + self.y_dimension_start, x + self.x_dimension_start)
                if self.sheet.cell_value(y + self.y_dimension_start, x + self.x_dimension_start) == '':
                    go_on = False
                    break
                else:
                    cell_string += \
                        f'{int(self.sheet.cell_value(y + self.y_dimension_start, x + self.x_dimension_start))}'
                    try:
                        cell_float.append(float(
                            f'{self.sheet.cell_value(y + self.y_dimension_start, x + self.x_dimension_start)}.0'))
                    except ValueError:
                        cell_float.append(float(
                            f'{self.sheet.cell_value(y + self.y_dimension_start, x + self.x_dimension_start)}'))
            if go_on:
                plaque_string_list.append(cell_string)
                plaque_float_list.append(cell_float)
        if self.is_simulation:
            parentals_string = []
            recombinants_string = []
            parentals_float = []
            recombinants_float = []
            for plaque in range(len(plaque_string_list)):
                plaque_string = plaque_string_list[plaque]
                plaque_float = plaque_float_list[plaque]
                if plaque_string == '11111111' or plaque_string == '00000000':
                    parentals_string.append(plaque_string)
                    parentals_float.append(plaque_float)
                else:
                    recombinants_string.append(plaque_string)
                    recombinants_float.append(plaque_float)
            segments_string = []
            segments_float = []
            for segment in range(8):
                segment_list_string = []
                segment_list_float = []
                for plaque in range(len(recombinants_string)):
                    segment_list_string.append(recombinants_string[plaque][segment])
                    segment_list_float.append(recombinants_float[plaque][segment])
                segments_string.append(segment_list_string)
                segments_float.append(segment_list_float)
            random_size = len(segments_string[0]) - 1
            generated_plaques_string = []
            generated_plaques_float = []
            for plaque in range(plaque_generation_size):
                if len(recombinants_string) == 0:
                    break
                random_list = []
                for segment in range(8):
                    try:
                        random_list.append(random.randint(0, random_size))
                    except ValueError:
                        random_list.append(0)
                new_plaque_string = ''
                new_plaque_float = []
                for segment in range(8):
                    random_value = random_list[segment]
                    new_plaque_string += segments_string[segment][random_value]
                    new_plaque_float.append(segments_float[segment][random_value])
                generated_plaques_string.append(new_plaque_string)
                generated_plaques_float.append(new_plaque_float)
            plaque_string_list = plaque_string_list + generated_plaques_string
            plaque_float_list = plaque_float_list + generated_plaques_float
        return plaque_string_list, plaque_float_list

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

    def get_recombinant_array(self, recombinants):
        segments = []
        for segment in range(8):
            segment_list = []
            for plaque in recombinants:
                segment_list.append(plaque[segment])
            segments.append(random.sample(segment_list, len(segment_list)))
        segments = np.array(segments).transpose()
        return segments

    def __str__(self):
        """Print dataframe of unique genotypes and appearances"""
        return f'{self.unique_genotypes_df}'


