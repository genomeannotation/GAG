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
    oldkeys = old.keys()
    newkeys = new.keys()
    for key in newkeys:
        if key not in oldkeys:
            return False
    return True

class StatsManager:

    increment_stats = ["seq_length", "num_genes", "num_mRNA", "num_exons", "num_introns", "num_CDS",\
            "CDS: complete", "CDS: start, no stop", "CDS: stop, no start", "CDS: no stop, no start",\
            "total_gene_length", "total_mRNA_length", "total_exon_length",\
            "total_intron_length", "total_CDS_length"]    
    min_stats = ["shortest_gene", "shortest_mRNA", "shortest_exon", "shortest_intron", "shortest_CDS"]
    max_stats = ["longest_gene", "longest_mRNA", "longest_exon", "longest_intron", "longest_CDS"]
    """
    calc_stats = ["mean gene length", "mean mRNA length", "mean exon length", "mean intron length",\
            "mean CDS length", "% of genome covered by genes", "% of genome covered by CDS",\
            "mRNAs per gene", "exons per mRNA", "introns per mRNA"]
            """
    calc_stats_formulae = {"mean gene length": ["total_gene_length", "num_genes"],\
            "mean mRNA length": ["total_mRNA_length", "num_mRNA"],\
            "mean exon length": ["total_exon_length", "num_exons"],\
            "mean intron length": ["total_intron_length", "num_introns"],\
            "mean CDS length": ["total_CDS_length", "num_CDS"]}
    calc_stats = ["mean gene length", "mean mRNA length", "mean exon length", "mean intron length",\
            "mean CDS length"]

    def __init__(self):
        self.ref_stats = {}
        self.alt_stats = {}
        self.initialize_dict(self.ref_stats)
        self.initialize_dict(self.alt_stats)

    def initialize_dict(self, d):
        for stat in self.increment_stats + self.min_stats + self.max_stats + self.calc_stats:
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

    def calculate_stat(self, stat):
        dividend_key = self.calc_stats_formulae[stat][0]
        divisor_key = self.calc_stats_formulae[stat][1]
        # Calculate for reference genome
        dividend = self.ref_stats[dividend_key]
        divisor = self.ref_stats[divisor_key]
        if divisor == 0:
            self.ref_stats[stat] = 0
        else:
            self.ref_stats[stat] = float(dividend) / float(divisor)

    def summary(self):
        for stat in self.calc_stats:
            self.calculate_stat(stat)
        stats_order = [key for keys in [self.increment_stats, self.min_stats, self.max_stats, self.calc_stats] for key in keys]
        return format_columns(["Reference Genome", "Modified Genome"], stats_order, [self.ref_stats, self.alt_stats], 5)



