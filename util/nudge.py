#!/usr/bin/env python
# coding=utf-8

import sys
import argparse


# Command line script to update .agp file
# According to a .bed file of locations
# which have been trimmed from the genome

def overlap(indices1, indices2):
    """Returns a boolean indicating whether two pairs of indices overlap."""
    if not (len(indices1) == 2 and len(indices2) == 2):
        return False
    if indices2[0] <= indices1[0] <= indices2[1]:
        return True
    elif indices2[0] <= indices1[1] <= indices2[1]:
        return True
    else:
        return False


def contains(indices1, indices2):
    """Returns a boolean indicating whether indices1 contain indices2"""
    if indices1[0] <= indices2[0] and indices1[1] >= indices2[1]:
        return True
    else:
        return False


def read_bed_file(filename):
    trimlist = []
    with open(filename, 'r') as bed:
        for line in bed:
            splitline = line.strip().split('\t')
            if len(splitline) != 3:
                return []
            else:
                try:
                    entry = [splitline[0], int(splitline[1]), int(splitline[2])]
                except ValueError:
                    sys.stderr.write("Error reading .bed file. Non-integer value ")
                    sys.stderr.write("in column 2 or 3. Here is the line:\n")
                    sys.stderr.write(line)
                    return []
                trimlist.append(entry)
        return trimlist


def fail_if_overlap(start, stop, trim_indices):
    if overlap([start, stop], trim_indices):
        sys.stderr.write("Collision! start/stop = %d/%d; \
                trim start/stop = %d/%d\n" %
                         (start, stop, trim_indices[0], trim_indices[1]))
        sys.exit()


def update_agp(agp_filename, trimlist):
    with open(agp_filename, 'r') as agp:
        for line in agp:
            fields = line.strip().split()
            seq_id = fields[0]
            start = int(fields[1])
            stop = int(fields[2])
            # In the case that there are multiple regions to trim in a single
            # sequence, trim from the end so indices don't get messed up
            to_trim_this_seq = [x for x in trimlist if x[0] == seq_id]
            to_trim_this_seq = sorted(to_trim_this_seq, key=lambda entry: entry[2], reverse=True)
            for trim_indices in to_trim_this_seq:
                trim_start = trim_indices[1]
                trim_stop = trim_indices[2]
                # Do nothing if trim indices fall after start/stop
                if trim_start > stop:
                    sys.stderr.write("trim region falls after stop; doing nothing\n")
                    continue
                # Fail spectacularly if trim region overlaps start->stop region
                fail_if_overlap(start, stop, trim_indices)
                # Fail spectacularly if trim region *contains* start->stop region
                if contains(trim_indices, [start, stop]):
                    sys.stderr.write("Collision! start/stop = %d/%d; \
                            trim start/stop = %d/%d\n" %
                                     (start, stop, trim_start, trim_stop))
                    sys.exit()
                # Trim region comes before our indices -- trim!
                bases_to_trim = trim_stop - trim_start + 1
                # Only adjust start if the trim region comes before it
                if trim_stop <= start:
                    start -= bases_to_trim
                stop -= bases_to_trim
                fields[1] = str(start)
                fields[2] = str(stop)
            print("\t".join(fields))


def update_gff(gff_filename, trimlist):
    with open(gff_filename, 'r') as gff:
        for line in gff:
            fields = line.strip().split()
            if len(fields) < 5:
                # A comment or something
                continue
            seq_id = fields[0]
            start = int(fields[3])
            stop = int(fields[4])
            # In the case that there are multiple regions to trim in a single
            # sequence, trim from the end so indices don't get messed up
            to_trim_this_seq = [x for x in trimlist if x[0] == seq_id]
            to_trim_this_seq = sorted(to_trim_this_seq, key=lambda entry: entry[2], reverse=True)
            for trim_indices in to_trim_this_seq:
                trim_start = trim_indices[1]
                trim_stop = trim_indices[2]
                # Do nothing if trim indices fall after start/stop
                if trim_start > stop:
                    sys.stderr.write("trim region falls after stop; doing nothing\n")
                    continue
                # Fail spectacularly if trim region overlaps start->stop region
                fail_if_overlap(start, stop, trim_indices)
                # Fail spectacularly if trim region *contains* start->stop region
                if contains(trim_indices, [start, stop]):
                    sys.stderr.write("Collision! start/stop = %d/%d; \
                            trim start/stop = %d/%d\n" %
                                     (start, stop, trim_start, trim_stop))
                    sys.exit()
                # Trim region comes before our indices -- trim!
                bases_to_trim = trim_stop - trim_start + 1
                # Only adjust start if the trim region comes before it
                if trim_stop <= start:
                    start -= bases_to_trim
                stop -= bases_to_trim
                fields[3] = str(start)
                fields[4] = str(stop)
            print("\t".join(fields))


def main():
    # Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed', required=True)
    parser.add_argument('-a', '--agp')
    parser.add_argument('-g', '--gff')
    args = parser.parse_args()
    # Read bed file (required)
    trimlist = read_bed_file(args.bed)
    if not trimlist:
        sys.stderr.write("Failed to read .bed file. Exiting...\n")
        sys.exit()
    # Read agp file, modify indices as needed
    if args.agp:
        update_agp(args.agp, trimlist)
    # Read gff file, modify indices as needed
    if args.gff:
        update_gff(args.gff, trimlist)


if __name__ == "__main__":
    main()
