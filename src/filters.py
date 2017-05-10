#!/usr/bin/env python
# coding=utf-8


# CDS

class MinCDSLengthFilter(object):
    def __init__(self, min_length=0):
        self.arg = min_length
        self.filter_mode = "REMOVE"
        return

    def apply(self, seq):
        count = 0
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.cds and mrna.cds.length() < self.arg:
                    mrna.death_flagged = True  # Destroy the mRNA that the cds lives on?
            to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
            count += len(to_remove)
            for mrna in to_remove:
                if self.filter_mode == "REMOVE":
                    print("Removing mRNA: " + mrna.identifier)
                    gene.remove_mrna(mrna.identifier)
                elif self.filter_mode == "FLAG":
                    print("Flagging mRNA: " + mrna.identifier)
                    mrna.cds.add_annotation('gag_flag', "cds_min_length:" + str(self.arg))
                elif self.filter_mode == "LIST":
                    print(mrna.identifier)
            for mrna in gene.mrnas:
                mrna.death_flagged = False
        if self.filter_mode == "REMOVE":
            print("\nRemoved " + str(count) + " mRNAs")
        elif self.filter_mode == "FLAG":
            print("\nFlagged " + str(count) + " mRNAs")
        elif self.filter_mode == "LIST":
            print(str(count) + " mRNAs")


class MaxCDSLengthFilter(object):
    def __init__(self, max_length=0):
        self.arg = max_length
        self.filter_mode = "REMOVE"
        return

    def apply(self, seq):
        count = 0
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.cds and 0 < self.arg < mrna.cds.length():
                    mrna.cds.add_annotation('gag_flag', "cds_max_length:" + str(self.arg))
                    mrna.death_flagged = True  # Destroy the mRNA that the cds lives on?
            to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
            count += len(to_remove)
            for mrna in to_remove:
                if self.filter_mode == "REMOVE":
                    print("Removing mRNA: " + mrna.identifier)
                    gene.remove_mrna(mrna.identifier)
                elif self.filter_mode == "FLAG":
                    print("Flagging mRNA: " + mrna.identifier)
                    mrna.cds.add_annotation('gag_flag', "cds_min_length:" + str(self.arg))
                elif self.filter_mode == "LIST":
                    print(mrna.identifier)
            for mrna in gene.mrnas:
                mrna.death_flagged = False
        if self.filter_mode == "REMOVE":
            print("\nRemoved " + str(count) + " mRNAs")
        elif self.filter_mode == "FLAG":
            print("\nFlagged " + str(count) + " mRNAs")
        elif self.filter_mode == "LIST":
            print(str(count) + " mRNAs")


# Exon

class MinExonLengthFilter(object):
    def __init__(self, min_length=0):
        self.arg = min_length
        self.filter_mode = "REMOVE"
        return

    def apply(self, seq):
        count = 0
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.exon and mrna.get_shortest_exon() < self.arg:
                    mrna.exon.add_annotation("gag_flag", "exon_min_length:" + str(self.arg))
                    mrna.death_flagged = True  # Destroy the mRNA that the exon lives on?
            to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
            count += len(to_remove)
            for mrna in to_remove:
                if self.filter_mode == "REMOVE":
                    print("Removing mRNA: " + mrna.identifier)
                    gene.remove_mrna(mrna.identifier)
                elif self.filter_mode == "FLAG":
                    print("Flagging mRNA: " + mrna.identifier)
                    mrna.cds.add_annotation('gag_flag', "cds_min_length:" + str(self.arg))
                elif self.filter_mode == "LIST":
                    print(mrna.identifier)
            for mrna in gene.mrnas:
                mrna.death_flagged = False
        if self.filter_mode == "REMOVE":
            print("\nRemoved " + str(count) + " mRNAs")
        elif self.filter_mode == "FLAG":
            print("\nFlagged " + str(count) + " mRNAs")
        elif self.filter_mode == "LIST":
            print(str(count) + " mRNAs")


class MaxExonLengthFilter(object):
    def __init__(self, max_length=0):
        self.arg = max_length
        self.filter_mode = "REMOVE"
        return

    def apply(self, seq):
        count = 0
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.exon and 0 < self.arg < mrna.get_longest_exon():
                    mrna.exon.add_annotation("gag_flag", "exon_max_length:" + str(self.arg))
                    mrna.death_flagged = True  # Destroy the mRNA that the exon lives on?
            to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
            count += len(to_remove)
            for mrna in to_remove:
                if self.filter_mode == "REMOVE":
                    print("Removing mRNA: " + mrna.identifier)
                    gene.remove_mrna(mrna.identifier)
                else:
                    print("Flagging mRNA: " + mrna.identifier)
                    mrna.cds.add_annotation('gag_flag', "cds_min_length:" + str(self.arg))
            for mrna in gene.mrnas:
                mrna.death_flagged = False
        if self.filter_mode == "REMOVE":
            print("\nRemoved " + str(count) + " mRNAs")
        else:
            print("\nFlagged " + str(count) + " mRNAs")


# Intron

class MinIntronLengthFilter(object):
    def __init__(self, min_length=0):
        self.arg = min_length
        self.filter_mode = "REMOVE"
        return

    def apply(self, seq):
        count = 0
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.exon and mrna.get_shortest_intron() < self.arg and mrna.get_shortest_intron() != 0:
                    mrna.exon.add_annotation("gag_flag", "intron_min_length:" + str(self.arg))
                    mrna.death_flagged = True  # Destroy the mRNA that the intron lives on?
            to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
            count += len(to_remove)
            for mrna in to_remove:
                if self.filter_mode == "REMOVE":
                    print("Removing mRNA: " + mrna.identifier)
                    gene.remove_mrna(mrna.identifier)
                elif self.filter_mode == "FLAG":
                    print("Flagging mRNA: " + mrna.identifier)
                    mrna.cds.add_annotation('gag_flag', "cds_min_length:" + str(self.arg))
                elif self.filter_mode == "LIST":
                    print(mrna.identifier)
            for mrna in gene.mrnas:
                mrna.death_flagged = False
        if self.filter_mode == "REMOVE":
            print("\nRemoved " + str(count) + " mRNAs")
        elif self.filter_mode == "FLAG":
            print("\nFlagged " + str(count) + " mRNAs")
        elif self.filter_mode == "LIST":
            print(str(count) + " mRNAs")


class MaxIntronLengthFilter(object):
    def __init__(self, max_length=0):
        self.arg = max_length
        self.filter_mode = "REMOVE"
        return

    def apply(self, seq):
        count = 0
        for gene in seq.genes:
            for mrna in gene.mrnas:
                if mrna.exon and 0 < self.arg < mrna.get_longest_intron():
                    mrna.exon.add_annotation("gag_flag", "intron_max_length:" + str(self.arg))
                    mrna.death_flagged = True  # Destroy the mRNA that the intron lives on?
            to_remove = [mrna for mrna in gene.mrnas if mrna.death_flagged]
            count += len(to_remove)
            for mrna in to_remove:
                if self.filter_mode == "REMOVE":
                    print("Removing mRNA: " + mrna.identifier)
                    gene.remove_mrna(mrna.identifier)
                elif self.filter_mode == "FLAG":
                    print("Flagging mRNA: " + mrna.identifier)
                    mrna.cds.add_annotation('gag_flag', "cds_min_length:" + str(self.arg))
                elif self.filter_mode == "LIST":
                    print(mrna.identifier)
            for mrna in gene.mrnas:
                mrna.death_flagged = False
        if self.filter_mode == "REMOVE":
            print("\nRemoved " + str(count) + " mRNAs")
        elif self.filter_mode == "FLAG":
            print("\nFlagged " + str(count) + " mRNAs")
        elif self.filter_mode == "LIST":
            print(str(count) + " mRNAs")


# Gene

class MinGeneLengthFilter(object):
    def __init__(self, min_length=0):
        self.arg = min_length
        self.filter_mode = "REMOVE"
        return

    def apply(self, seq):
        count = 0
        for gene in seq.genes:
            if gene.length() < self.arg:
                gene.death_flagged = True  # Destroy the gene?
        for gene in seq.genes:
            if gene.death_flagged:
                count += 1
                if self.filter_mode == "REMOVE":
                    print("Removing gene: " + gene.identifier)
                elif self.filter_mode == "FLAG":
                    print("Flagging gene: " + gene.identifier)
                    gene.add_annotation("gag_flag", "gene_min_length:" + str(self.arg))
                elif self.filter_mode == "LIST":
                    print(gene.identifier)
        if self.filter_mode == "REMOVE":
            seq.genes = [gene for gene in seq.genes if not gene.death_flagged]
        for gene in seq.genes:
            gene.death_flagged = False
        if self.filter_mode == "REMOVE":
            print("\nRemoved " + str(count) + " genes")
        elif self.filter_mode == "FLAG":
            print("\nFlagged " + str(count) + " genes")
        elif self.filter_mode == "LIST":
            print(str(count) + " genes")


class MaxGeneLengthFilter(object):
    def __init__(self, max_length=0):
        self.arg = max_length
        self.filter_mode = "REMOVE"
        return

    def apply(self, seq):
        count = 0
        for gene in seq.genes:
            if 0 < self.arg < gene.length():
                gene.add_annotation("gag_flag", "gene_max_length:" + str(self.arg))
                gene.death_flagged = True  # Destroy the gene?
        for gene in seq.genes:
            if gene.death_flagged:
                count += 1
                if self.filter_mode == "REMOVE":
                    print("Removing gene: " + gene.identifier)
                elif self.filter_mode == "FLAG":
                    print("Flagging gene: " + gene.identifier)
                    gene.add_annotation("gag_flag", "gene_min_length:" + str(self.arg))
                elif self.filter_mode == "LIST":
                    print(gene.identifier)
        if self.filter_mode == "REMOVE":
            seq.genes = [gene for gene in seq.genes if not gene.death_flagged]
        for gene in seq.genes:
            gene.death_flagged = False
        if self.filter_mode == "REMOVE":
            print("\nRemoved " + str(count) + " genes")
        elif self.filter_mode == "FLAG":
            print("\nFlagged " + str(count) + " genes")
        elif self.filter_mode == "LIST":
            print(str(count) + " genes")
