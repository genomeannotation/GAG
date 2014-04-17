#!/usr/bin/env python


class SeqFixer:

    def __init__(self):
        self.terminal_ns = False
        self.start_stop_codons = False
        self.dirty = False

    def fix_terminal_ns(self):
        self.terminal_ns = True
        self.dirty = True

    def fix_start_stop_codons(self):
        self.start_stop_codons = True
        self.dirty = True

    def fix(self, seq):
        if self.terminal_ns:
            seq.remove_terminal_ns()
        if self.start_stop_codons:
            seq.create_starts_and_stops()

