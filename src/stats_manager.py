#!/usr/bin/env python

def format_column(column, spacing):
    # First, get the uniform length
    longest = 0
    for item in column:
        length = len(item)+spacing
        if length > longest:
            longest = length
    # Now format
    return [item+(' '*(longest-len(item))) for item in column]

def format_columns(column_names, key_order, dicts, spacing = 3):
    # Build key column
    columns = [['', '']]
    columns[0].extend(key_order)
    columns[0] = format_column(columns[0], spacing)
    
    # Notes: Python automatically sorts dictionary contents by key, so this will work.
    # TODO: Nevertheless, make this code less hacky
    
    for i, dic in enumerate(dicts):
        new_column = [column_names[i], '-'*len(column_names[i])]
        new_column.extend([str(dic[key]) for key in key_order])
        columns.append(format_column(new_column, spacing))
        
    # Finally, stringify the table
    tbl_str = ''
    for i in range(len(columns[0])): # For each row
        for column in columns: # For each column
            tbl_str += column[i]
        tbl_str += '\n'
    return tbl_str

def validate_dicts(old, new):
    if len(old) != len(new):
        return False
    oldkeys = old.keys()
    newkeys = new.keys()
    for key in oldkeys:
        if key not in newkeys:
            return False
    return True

class StatsManager:

    increment_stats = ["seq_length", "num_genes", "num_mRNA", "num_exons", "num_CDS",\
            "CDS: complete", "CDS: start, no stop", "CDS: stop, no start", "CDS: no stop, no start",\
            "total_gene_length", "total_mRNA_length", "total_exon_length",\
            "total_intron_length", "total_CDS_length"]    
    min_stats = ["shortest_gene", "shortest_mRNA", "shortest_exon", "shortest_intron", "shortest_CDS"]
    max_stats = ["longest_gene", "longest_mRNA", "longest_exon", "longest_intron", "longest_CDS"]

    def __init__(self):
        self.ref_stats = {}
        self.alt_stats = {}
        self.initialize_dict(self.ref_stats)
        self.initialize_dict(self.alt_stats)

    def initialize_dict(self, d):
        for stat in self.increment_stats + self.min_stats + self.max_stats:
            d[stat] = 0

    def clear_alt(self):
        self.initialize_dict(self.alt_stats)

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
        stats_order = [key for keys in [self.increment_stats, self.min_stats, self.max_stats] for key in keys]
        return format_columns(["Reference Genome", "Modified Genome"], stats_order, [self.ref_stats, self.alt_stats], 5)



