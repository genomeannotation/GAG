#!/usr/bin/env python

from src.annotator import Annotator

class Genome:

    """Genome is the master class for storing all of the data in a genome.
    It owns and provides access to the fasta and genes.
    
    """

    def __init__(self):
        self.fasta = None
        self.genes = []
        self.annot = None
        self.entries = []

    # this also removes empty genes; could use a better name maybe...
    def remove_mrnas_with_cds_shorter_than(self, min_length):
        """Removes all mRNAs containing CDS shorter than min_length. If this causes
        a gene to contain no mRNAs, the gene will also be removed.

        """
        if self.genes:
            to_remove = []
            for gene in self.genes:
                gene.remove_mrnas_with_cds_shorter_than(min_length)
                if not gene.mrnas:
                    to_remove.append(gene)
            for g in to_remove:
                self.genes.remove(g)

    def remove_first_cds_segment_if_shorter_than(self, min_length):
        """Removes only the first CDS segment if it's shorter than min_length. If
        this causes a gene to contain no mRNAs, the gene will also be removed.

        """
        if self.genes:
            for gene in self.genes:
                gene.remove_first_cds_segment_if_shorter_than(min_length)

    # maybe should be called 
    #  'create_start_codon_GenePart if sequence contains start codon'
    def create_starts_and_stops(self):
        """Checks the actual sequence to see if there are start and stop codons
        on each gene. If there are, the appropriate GeneParts are made and added to
        the gene.

        """

        for gene in self.genes:
            gene.create_starts_and_stops(self.fasta)

    def write_string(self, genes = None, errors = None):
        """Writes this genome to NCBI's tbl format and returns it as a string.

        """

        output = ''

        if self.fasta == None or self.genes == None:
            return output
        if self.annot == None:
            self.annot = Annotator()
       
        output += '>Feature SeqId\n'

        for seq in self.fasta.entries:
            if genes != None and not genes:
                return output

            entries = []
            for gene in self.genes:
                if gene.seq_name != seq[0]:
                    continue

                if genes != None and gene.name not in genes:
                    continue

                if genes != None:
                    genes.remove(gene.name)

                newEntries = gene.to_tbl_entries(self.annot)
                for entry in newEntries:
                    entries.append(entry)

            # If there are any entries, write this section of the tbl file
            if len(entries) > 0:
                output += '>Feature '+seq[0]+'\n'
                output += '1\t'+str(len(seq[1]))+\
                          '\tREFERENCE\n\t\t\tPBARC\t12345\n'

                for entry in entries:    
                    output += entry.write_to_string()+'\n'
        return output

    def remove_all_gene_segments(self, prefix):
        """Remove all genes whose names contain the prefix string.

        Note: Prefix is not actually a prefix. If the gene name contains the string
        at all, it will be removed.

        """
        if len(prefix) > 0:
            self.genes = [g for g in self.genes if not prefix in g.identifier]

    def remove_genes_marked_for_removal(self):
        """Removes all genes previously marked for removal, by setting their indices at
        [0, 0]

        """
        for gene in reversed(self.genes):
            if gene.indices[0] == 0 and gene.indices[1] == 0:
                self.genes.remove(gene)

    def get_locus_tag(self):
        firstgene = self.genes[0]
        gene_id = str(firstgene.identifier)
        locus_tag = gene_id.split('_')[0]
        return locus_tag

    def rename_maker_mrnas(self):
        """Renames weird maker mRNA ids to ids in the format <locus_tag>1000000

        """

        locus_tag = self.get_locus_tag()
        count = 1000000
        for gene in self.genes:
            for mrna in gene.mrnas:
                if mrna.is_maker_mrna():
                    old_name = mrna.identifier
                    new_name = locus_tag + '_' + str(count)
                    mrna.identifier = new_name
                    self.annot.rename_mrna(old_name, new_name)
                    count += 1

    def invalidate_region(self, seq, start, stop):
        """Invalidates a region on a specified sequence. Any CDS or Exons containing the
        region will be trimmed, and any CDS or Exons encompassing the region will be
        removed.

        """
        
        for gene in self.genes:
            if gene.seq_name == seq:
                gene.invalidate_region(start, stop)

    def trim_region(self, seq, start, stop):
        """Removes a region from a specified sequence. The GFF indices will be adjusted
        accordingly.

        """

        if self.fasta:
            self.fasta.trim_seq(seq, start, stop)
        if self.genes:
            offset = -(stop - start + 1)
            for gene in self.genes:
                if gene.seq_name == seq:
                    gene.adjust_indices(offset, stop)

    def remove_seq(self, seq_id, force=False):
        """Removes a specified sequence. Fails gracefully if the sequence contains
        features. Can be forced to remove the sequence regardless.

        """

        if force or not self.contains_gene_on_seq(seq_id):
            self.fasta.remove_seq(seq_id)
        else:
            print("Sorry, that sequence contains features. \
                    Try '-F' to force removal\n")

    def get_gene_seq_info(self, gene_name):
        """Fetch the name and indices of a gene as a tuple.

        """

        for gene in self.genes:
            if gene.name == gene_name:
                return [gene.seq_name, gene.indices]

    def check_gene_for_invalid_begin_or_end(self, gene_name):
        """Checks a gene for Ns at the beginning or end. If there are Ns, the gene
        indices will be trimmed to not contain them.

        """

        seqinfo = self.get_gene_seq_info(gene_name)
        if not seqinfo:
            return
        seq_id = seqinfo[0]
        begin = seqinfo[1][0]
        end = seqinfo[1][1]
        bases_forward = self.fasta.how_many_Ns_forward(seq_id, begin)
        if bases_forward != 0:
            self.invalidate_region(seq_id, begin, begin+bases_forward-1)
        bases_backward = self.fasta.how_many_Ns_backward(seq_id, end)
        if bases_backward != 0:
            self.invalidate_region(seq_id, end-bases_backward+1, end)

    def contains_gene_on_seq(self, seq_id):
        for gene in self.genes:
            if gene.seq_name == seq_id:
                return True
        return False
