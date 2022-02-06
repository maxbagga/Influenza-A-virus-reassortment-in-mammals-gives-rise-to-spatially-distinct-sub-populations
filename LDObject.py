"""Defines object used in linkage disequilibrium calculations"""

class LDObject:
    def __init__(self, segment_values):
        """Takes in a list of values for one segment"""
        self.segment_values = segment_values
        self.prob1 = self.probability_virus(1.0)
        self.prob2 = self.probability_virus(0.0)
        self.hprob1 = self.prob1**2
        self.hprob2 = self.prob2**2
        self.outcome1 = self.prob1*len(segment_values)
        self.outcome2 = self.prob2*len(segment_values)

    def probability_virus(self, value):
        """Calculates probability of getting a certain parent for segment"""
        count = 0
        for plaque in self.segment_values:
            if plaque == value:
                count += 1
        return count/len(self.segment_values)
