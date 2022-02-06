from scipy.stats import f_oneway


def dose_ANOVA(low_dose, high_dose):
    high_dose_days_swi = []
    low_dose_days_swi = []
    for individual in low_dose.individuals:
        for day in individual.days:
            low_dose_days_swi.append(day.genotype_analyzer.shannon_wiener)
    for individual in high_dose.individuals:
        for day in individual.days:
            high_dose_days_swi.append(day.genotype_analyzer.shannon_wiener)
    print(f_oneway(low_dose_days_swi, high_dose_days_swi))