#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import cmd
import readline
import sys
from src.console_controller import ConsoleController


class GagCmd(cmd.Cmd):

## Setup, loading and saving sessions, exit, etc.

    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = "GAG> "
        readline.set_history_length(1000)
        readline.read_history_file('.gaghistory')
        self.controller = ConsoleController() 

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def postcmd(self, stop, line):
        if hasattr(self, 'output') and self.output:
            print self.output
            self.output = None
        return stop

    def parseline(self, line):
        if '|' in line:
            commands = ['']
            quotes = ''
            for c in line:
                if c == '|' and quotes == '':
                    commands[-1] = commands[-1].strip()
                    commands.append('')
                else:
                    if quotes == '' and (c == "'" or c == '"'):
                        quotes = c
                    elif quotes == c:
                        quotes = ''
                    commands[-1] += c
            commands[-1] = commands[-1].strip()
            if len(commands) > 1:
                return 'pipe', commands, line
        return cmd.Cmd.parseline(self, line)

    def do_pipe(self, args):
        buf = None
        for arg in args:
            if buf:
                self.controller.input = buf
            else:
                self.controller.input = ''
            self.onecmd(arg)
            if hasattr(self, 'output') and self.output:
                buf = self.output

    def do_barfsession(self, line):
        self.controller.barf_session(line)

    def do_loadsession(self, line):
        self.controller.load_session(line)

    def help_exit(self):
        print("Exit this console.\n")

    def do_exit(self, line):
        return True

    def do_ls(self, line):
        self.output = self.controller.ls(line)

    def do_cat(self, line):
        self.output = self.controller.cat(line)

    def do_grep(self, line):
        self.output = self.controller.grep(line)

    def do_sed(self, line):
        self.output = self.controller.sed(line)


## Reading in files

    def help_readfasta(self):
        print("Usage: readfasta <file_name>\n")
        print("Read the fasta file. Any unsaved changes")
        print("to the currently loaded fasta will be lost.\n")

    def do_readfasta(self, line):
        self.controller.read_fasta(line)

    def help_readgff(self):
        print("Usage: readgff <file_name>\n")
        print("Read the gff file. Any unsaved changes")
        print("to the currently loaded gff will be lost.\n")

    def do_readgff(self, line):
        self.controller.read_gff(line)

    def help_readtrinotate(self):
        print("Usage: readtrinotate <file_name>\n")
        print("Read annotations from a trinotate file.\n")

    def do_readtrinotate(self, line):
        self.controller.read_trinotate(line)


## Manipulate genome

    def help_applybed(self):
        print("Usage: applybed <file_name>\n")
        print("Applies a bed file to the data. This will")
        print("trim out any sequences that aren't in the")
        print("ranges of the bed file and automagically")
        print("update the GFF file accordingly.\n")

    def do_applybed(self, line):
        self.controller.apply_bed(line)

    def do_ducttapeseqframes(self, line):
        self.controller.duct_tape_seq_frames(line)


## Output info to console

    def do_barfgenegff(self, line):
        self.output = self.controller.barf_gff(line)

    def do_barfseq(self, line):
        self.output = self.controller.barf_seq(line)

    def do_barfgenetbl(self, line):
        self.output = self.controller.barf_gene_tbl(line)


## Output info to file

    def help_writetbl(self):
        print("Usage: writetbl <file_name>\n")
        print("Write a sweet feature table to the specified file.\n")

    def do_writetbl(self, line):
        self.controller.write_tbl(line)


########################################################################

if __name__ == '__main__':
    GagCmd().cmdloop("Welcome to the GAG console!")
