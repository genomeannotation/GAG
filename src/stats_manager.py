#!/usr/bin/env python

def validate_dicts(old, new):
    # TODO check keys?
    return len(old) == len(new)

class StatsManager:

    increment_stats = ["seq_length", "num_genes", "num_mRNA", "num_CDS",\
            "total_gene_length", "total_mRNA_length", "total_CDS_length"]    
    min_stats = ["shortest_gene", "shortest_mRNA", "shortest_CDS"]
    max_stats = ["longest_gene", "longest_mRNA", "longest_CDS"]

    def __init__(self):
        self.ref_stats = {}
        self.alt_stats = {}
        self.initialize_dict(self.ref_stats)
        self.initialize_dict(self.alt_stats)

    def initialize_dict(self, d):
        for stat in self.increment_stats + self.min_stats + self.max_stats:
            d[stat] = 0

    def clear_ref(self):
        self.initialize_dict(self.ref_stats)

    def update_ref(self, stats):
        self.update_stats(self.ref_stats, stats)

    def update_alt(self, stats):
        self.update_stats(self.alt_stats, stats)

    def update_stats(self, old, new):
        if not validate_dicts(old, new):
            return
        for stat in self.increment_stats:
            old[stat] += new[stat]
        for stat in self.min_stats:
            if old[stat] == 0:
                old[stat] = new[stat]
            elif new[stat] < old[stat] and new[stat] != 0:
                old[stat] = new[stat]
        for stat in self.max_stats:
            if new[stat] > old[stat]:
                old[stat] = new[stat]

    def summary(self):
        result = "\t\tReference Genome\tModified Genome\n"
        result += "\t\t----------------\t---------------\n"
        for stat in self.increment_stats:
            # Subtract a tab for longer stat names #KindaGhetto
            if len(stat) > 10:
                result += stat + ":\t"
            else: result += stat + ":\t\t"
            result += str(int(self.ref_stats[stat]))
            result += "\t\t" + str(int(self.alt_stats[stat])) + "\n"
        for stat in self.min_stats:
            result += stat + ":\t\t" + str(int(self.ref_stats[stat]))
            result += "\t\t" + str(int(self.alt_stats[stat])) + "\n"
        for stat in self.max_stats:
            result += stat + ":\t\t" + str(int(self.ref_stats[stat]))
            result += "\t\t" + str(int(self.alt_stats[stat])) + "\n"
        return result



