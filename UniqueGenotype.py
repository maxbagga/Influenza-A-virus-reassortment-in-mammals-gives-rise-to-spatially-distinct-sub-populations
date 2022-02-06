from scipy.stats import binom as b
import numpy as np
import pandas as pd
import statistics
"""Genotype object that goes in genotype analyzer"""


class Genotype:
    def __init__(self, genotype, appearances, ldos, number_of_plaques):
        """Genotype object that goes in genotype analyzer"""
        self.genotype = genotype
        self.appearances = appearances
        self.ldos = ldos
        self.number_of_plaques = number_of_plaques
        self.prob, self.expected, self.id = self.get_info()
        self.prob_of_app = b.pmf(self.appearances, self.number_of_plaques, self.prob)
        self.pval, self.is_less = self.get_p()

    def get_info(self):
        prob = 1
        id = ""
        for segment in range(len(self.genotype)):
            if self.genotype[segment] == 1:
                prob *= self.ldos[segment].prob1
            else:
                prob *= self.ldos[segment].prob2
            id += f'{self.genotype[segment]}'
        return prob, prob * self.number_of_plaques, id

    def get_p(self):

        if self.appearances < self.expected:
            # return 2 * b.cdf(self.appearances, self.number_of_plaques, self.prob), True
            pval = 2*b.cdf(self.appearances, self.number_of_plaques, self.prob), True
        else:
            # return 2 * (1 - b.cdf(self.appearances, self.number_of_plaques, self.prob)), False
            pval = 2*(1-b.cdf(self.appearances-1, self.number_of_plaques, self.prob)), False
        if pval[0] > 1.0:
            return 1.0, pval[1]
        else:
            return pval

    def __str__(self):
        return \
            f'genotype: {self.id} Appearances: {self.appearances}  Expected Appearances: {round(self.expected, 5)} ' \
                f'Probability: {round(self.prob, 5)} Probability of Appearances: {round(5, self.prob_of_app)}'


class GenotypeAnalyzer:
    def __init__(self, Genotypes):
        """Genotype analyzer object that is used to obtain diversity statistics"""
        self.Genotypes = Genotypes
        self.pn1ss = self.get_pnss_wo_parentals(Genotypes)
        self.pn2ss = self.get_pnss_wo_parentals_other(Genotypes)
        self.differences, self.appearances, self.expecteds, self.ids, self.probabilities, self.prob_of_app, \
            self.shannon_wiener, self.pvals, self.is_lesser, self.ldoss, self.number_of_plaquess, self.unique, \
            self.shannon_wiener_wop, self.sum_p, self.unique_wop, self.shannon_wiener_wop_woo, self.sum_p1, \
            self.unique_wop_woo, self.plaques, self.plaques_wop, self.plaques_wop_woo = self.get_info()
        self.avg_diff = statistics.mean(self.differences)
        self.stdev_diff = statistics.pstdev(self.differences)
        self.df = self.get_df()
        self.avg_poa = statistics.mean(self.prob_of_app)
        try:
            self.stdev_poa = statistics.pstdev(self.prob_of_app)
        except AssertionError:
            self.stdev_poa = 0
        self.max_diversity = np.log(self.unique)
        self.max_diversity_wop = np.log(self.unique_wop)
        self.max_diversity_wop_woo = np.log(self.unique_wop_woo)
        self.evenness = self.shannon_wiener/self.max_diversity
        self.evenness_wop = self.shannon_wiener_wop/self.max_diversity_wop
        self.evenness_wop_woo = self.shannon_wiener_wop_woo/self.max_diversity_wop_woo
        self.avg_pval = statistics.mean(self.pvals)
        self.stdev_pval = statistics.pstdev(self.pvals)
        self.df_stats = self.get_df_stats()

    def get_info(self):
        differences = []
        appearances = []
        expecteds = []
        ids = []
        probabilities = []
        prob_of_apps = []
        shannon_wiener = 0
        shannon_wiener_wop = 0
        shannon_wiener_wop_woo = 0
        sum_p = 0
        sum_p1 = 0
        pvals = []
        is_lesser = []
        ldoss = []
        number_of_plaquess = []
        count = 0
        # for Genotype in self.Genotypes:
        #     if Genotype.appearances != 0:
        #         count += 1
        count_wop = 0
        count_wop_woo = 0
        plaques = 0
        plaques_wop = 0
        plaques_wop_woo = 0
        for genotype in self.Genotypes:
            differences.append(genotype.appearances - genotype.expected)
            appearances.append(genotype.appearances)
            expecteds.append(genotype.expected)
            ids.append(genotype.id)
            probabilities.append(genotype.prob)
            prob_of_apps.append(genotype.prob_of_app)
            number_of_plaquess.append(genotype.number_of_plaques)
            pn = genotype.appearances / number_of_plaquess[0]
            try:
                pn1 = genotype.appearances / self.pn1ss
            except ZeroDivisionError:
                pn1 = 0
            try:
                pn2 = genotype.appearances / self.pn2ss
            except ZeroDivisionError:
                pn2 = np.NaN
            if genotype.appearances != 0:
                shannon_wiener += pn*np.log(pn)
                count += 1
                plaques += genotype.appearances
            if genotype.appearances != 0 and genotype.id != "11111111" and genotype.id != "00000000":
                shannon_wiener_wop += pn1 * np.log(pn1)
                sum_p += pn1
                count_wop += 1
                plaques_wop += genotype.appearances
            if genotype.appearances != 0 and genotype.id != "11111111" and genotype.id != "00000000" and \
                    genotype.id != "11110111":
                shannon_wiener_wop_woo += pn2 * np.log(pn2)
                sum_p1 += pn2
                count_wop_woo += 1
                plaques_wop_woo += genotype.appearances
            pvals.append(genotype.pval)
            is_lesser.append(genotype.is_less)
            ldoss.append(genotype.ldos)
        return differences, appearances, expecteds, ids, probabilities, prob_of_apps, -1 * shannon_wiener, pvals, \
               is_lesser, ldoss, number_of_plaquess, count, -1 * shannon_wiener_wop, sum_p, count_wop, -1 * \
               shannon_wiener_wop_woo, sum_p1, count_wop_woo, plaques, plaques_wop, plaques_wop_woo

    def get_df(self):
        df = pd.DataFrame()
        df['Genotype:'] = self.ids
        df['Appearances:'] = self.appearances
        df['Expected Appearances:'] = self.expecteds
        df['Difference:'] = self.differences
        df['Probability:'] = self.probabilities
        df['Probability of Appearances:'] = self.prob_of_app
        df['p-Value:'] = self.pvals
        return df

    def get_df_stats(self):
        df = pd.DataFrame()
        df['Unique Genotypes:'] = [self.unique]
        df['Average Difference:'] = [self.avg_diff]
        df['STDEV of Differences:'] = [self.stdev_diff]
        df['Average Probability of Appearances:'] = [self.avg_poa]
        df['STDEV of Probability of Appearances'] = [self.stdev_poa]
        df['Shannon-Wiener:'] = [self.shannon_wiener]
        df['Maximum Diversity:'] = [self.max_diversity]
        df['Evenness:'] = [self.evenness]
        df['Average p-Value:'] = [self.avg_pval]
        df['STDEV of p-Values:'] = [self.stdev_pval]
        return df

    def get_pvals(self):
        return self.ids, self.pvals, self.is_lesser

    def get_pnss_wo_parentals(self, genotypes):
        pn1ss = genotypes[0].number_of_plaques
        for genotype in genotypes:
            if genotype.id == "11111111" or genotype.id == "00000000":
                pn1ss -= genotype.appearances
        return pn1ss

    def get_pnss_wo_parentals_other(self, genotypes):
        pn1ss = genotypes[0].number_of_plaques
        for genotype in genotypes:
            if genotype.id == "11111111" or genotype.id == "00000000" or genotype.id == "11110111":
                pn1ss -= genotype.appearances
        return pn1ss

    def __str__(self):
        return f'Unique Genotypes = {self.unique} Average Difference: {round(self.avg_diff, 5)} ' \
            f'STDEV of Differences: {round(self.stdev_diff, 5)} Average Probability of Appearances: ' \
            f'{round(self.avg_poa, 5)} STDEV of Probability of Appearances {round(self.stdev_poa, 5)} ' \
            f'Shannon-Wiener: {round(self.shannon_wiener, 5)} Maximum Diversity: {round(self.max_diversity, 5)} ' \
            f'Evenness: {round(self.evenness, 5)} Average p-Value: {round(self.avg_pval, 5)} STDEV of p-Values ' \
            f'{round(self.stdev_pval, 5)}'

    def __add__(self, other):
        cumulative = [[[], []], [[], []]]
        temp = other.Genotypes[:]
        for genotype1 in range(len(self.Genotypes)):
            in_both = False
            for genotype2 in range(len(other.Genotypes)):
                if self.Genotypes[genotype1].genotype == other.Genotypes[genotype2].genotype:
                    temp.remove(other.Genotypes[genotype2])
                    cumulative[0][0].append(self.Genotypes[genotype1].genotype)
                    cumulative[0][1].append(self.appearances[genotype1])
                    cumulative[1][0].append(other.Genotypes[genotype2].genotype)
                    cumulative[1][1].append(other.appearances[genotype2])
                    in_both = True
            if not in_both:
                cumulative[0][0].append(self.Genotypes[genotype1].genotype)
                cumulative[0][1].append(self.appearances[genotype1])
                cumulative[1][0].append(self.Genotypes[genotype1].genotype)
                cumulative[1][1].append(0)
        for genotype1 in temp:
            for genotype2 in range(len(other.Genotypes)):
                if genotype1.genotype == other.Genotypes[genotype2].genotype:
                    cumulative[0][0].append(other.Genotypes[genotype2].genotype)
                    cumulative[0][1].append(0)
                    cumulative[1][0].append(other.Genotypes[genotype2].genotype)
                    cumulative[1][1].append(other.appearances[genotype2])
        Genotypes1 = []
        Genotypes2 = []
        for genotype in range(len(cumulative[0][0])):
            Genotypes1.append(Genotype(cumulative[0][0][genotype], cumulative[0][1][genotype], self.ldoss[0],
                                       self.number_of_plaquess[0]))
            Genotypes2.append(Genotype(cumulative[1][0][genotype], cumulative[1][1][genotype], other.ldoss[0],
                                       other.number_of_plaquess[0]))
        return GenotypeAnalyzer(Genotypes1), GenotypeAnalyzer(Genotypes2)

