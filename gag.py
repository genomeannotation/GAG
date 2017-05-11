#!/usr/bin/env python
# coding=utf-8
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import argparse

from src.controller import Controller


def main():
    version = "v1.1.1"
    parser = argparse.ArgumentParser(
        epilog="""
    Docs at http://genomeannotation.github.io/GAG/
    Bugs and feature requests at https://github.com/genomeannotation/GAG/issues
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-f', '--fasta', required=True)
    parser.add_argument('-g', '--gff', required=True)
    parser.add_argument('-v', '--version', action='version', version="GAG " + version)
    parser.add_argument('-a', '--anno')
    parser.add_argument('-t', '--trim')
    parser.add_argument('-o', '--out')
    parser.add_argument('--fix_start_stop', action='store_true')
    parser.add_argument('--fix_terminal_ns', action='store_true')
    parser.add_argument('-rcs', '--remove_cds_shorter_than')
    parser.add_argument('-rcl', '--remove_cds_longer_than')
    parser.add_argument('-res', '--remove_exons_shorter_than')
    parser.add_argument('-rel', '--remove_exons_longer_than')
    parser.add_argument('-ris', '--remove_introns_shorter_than')
    parser.add_argument('-ril', '--remove_introns_longer_than')
    parser.add_argument('-rgs', '--remove_genes_shorter_than')
    parser.add_argument('-rgl', '--remove_genes_longer_than')
    parser.add_argument('-fcs', '--flag_cds_shorter_than')
    parser.add_argument('-fcl', '--flag_cds_longer_than')
    parser.add_argument('-fes', '--flag_exons_shorter_than')
    parser.add_argument('-fel', '--flag_exons_longer_than')
    parser.add_argument('-fis', '--flag_introns_shorter_than')
    parser.add_argument('-fil', '--flag_introns_longer_than')
    parser.add_argument('-fgs', '--flag_genes_shorter_than')
    parser.add_argument('-fgl', '--flag_genes_longer_than')
    parser.add_argument('-ses', '--skip_empty_scaffolds', action='store_true')
    parser.add_argument('-gct', '--genome_center_tag', default='ncbi')
    parser.add_argument('-wpf', '--wgs_accession_prefix')
    parser.add_argument('-ref', '--reference_qualifier', default='PBARC:12345')
    parser.add_argument('-tif', '--transcript_id_format', default='suffix', choices=['prefix', 'suffix'])
    args = parser.parse_args()
    controller = Controller()
    controller.execute(args)


if __name__ == '__main__':
    main()
