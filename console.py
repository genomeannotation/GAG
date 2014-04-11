#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import cmd
import readline
import sys
import traceback
from src.console_controller import ConsoleController

def try_catch(command, args):
    try:
        if args is not None:
            return command(*args)
        else:
            return command()
    except:
        print("Sorry, that command raised an exception. Here's what I know:\n")
        print(traceback.format_exc())


class GagCmd(cmd.Cmd):

## Setup, loading and saving sessions, exit, etc.

    intro = "................................................................   .       .    \n........................................,,........................  .... .  ....\n..............~++=:,............... ..=OMNZ:....... ....... ..:=?77I=,... .. ...\n...........:7DMMMMMMD7:............. ~DMMMMZ...............=$DMMMMMMMM8+........\n.........:7NMMMMMMMMMMN7,............$MMMMMN=. ....... ..~ZMMMMMMMMMMMMZ,.   ...\n........?NMMMNZ+~+IONMMN?....... ...INMMNMMMI .. .......INMMMDZ?~,,:+I?:. . .  .\n.......+NMMN7:  ... :?$I: .........INMMNZMMMZ.........:$MMMD?........   ........\n......:8MMM?..........  ..........7MMMD+=NMM8,... .. :OMMMO~ ...................\n......7MMMZ,.....................7MMMD= ,8MMM? ......ZMMMZ,.....................\n.....:8MMN+....................,$MMMN+.. IMMMO......~NMMD~ ....,+$ZI+==~,.......\n.....~DMMD~....~7OZ$Z8NNZ:....+DMMMD$I?~:~DMMN+. ...+NMM8......IMMMMMMMMO: ..  .\n......$MMM8:. ,OMMMMMMMMN=...INMMN$$MMMMMMMMMMZ.....~DMMMI.....~8MMMMMMMZ: ...  \n......,ZMMMD?,.?NMMMMMNZ~. .INMMD+.=ONMMMMMMMMN~. ...+DMMM$:.  ..,=$MMMD: ...   \n.......,$MMMMM87+IONMMMI...?NMMN= ....:=+??OMMM$..... =8MMMM87=,.,+8MMMZ.  .   .\n.........+8MMMMMMMMMMMO~. ~DMMN?...... . ..=NMM$.......:$NMMMMMMNMMMMMO:.... . .\n...........:I8NMMMMNO+,...:OMD?. ....... ...+Z7~.........:I8NMMMMMMMDI,    .   .\n............. .:~~:. ......,,,.. . ... . . . ..  . ... . .  .~?7$7?~.        .  \n................................ . . .   . . .   . . .   . . .         .       .\n................. ....... ......   . .     . .     . .     . .         .       ."+\
    "\nWelcome to the GAG console!\n"+\
    "Type 'help' for available commands.\n"

    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = "GAG> "
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")
        self.controller = ConsoleController() 

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_barffolder(self):
        print("Usage: barffolder <directory>\n")
        print("Writes gff, fasta and trinotate files to the specified directory.\n")

    def do_barffolder(self, line):
        print(try_catch(self.controller.barf_folder, [line]))

    def help_load(self):
        print("This command takes you the GAG LOAD menu. There you can specify the location of")
        print("your files and load them into memory.")
        print("Alternately, just type 'load <path>' and avoid the submenu altogether.")

    def do_load(self, line):
        path_to_load = line.strip()
        loadcmd = LoadCmd(self.prompt, self.controller, path_to_load)
        loadcmd.cmdloop()

    def help_exit(self):
        print("Exit this console.\n")

    def do_exit(self, line):
        return True
        
    def help_setfilterarg(self):
        print("\nUsage: modifyfilterarg <filter_name> <filter_arg> <value>\n\nSets a specified filter argument to a specified value.\n")
    
    def do_setfilterarg(self, line):
        try_catch(self.controller.set_filter_arg, line.split(' '))
    
    def help_getfilterarg(self):
        print("\nUsage: getfilterarg <filter_name> <filter_arg> <value>\n\nOutputs the value of the specified filter argument.\n")
    
    def do_getfilterarg(self, line):
        print(try_catch(self.controller.get_filter_arg, line.split(' ')))
        
    def help_getfilterhelp(self):
        print("\nUsage: getfilterhelp [filter_name]\n\nGets help for a filter. Omit filter name to list info for all filters.\n")
        
    def do_getfilterhelp(self, line):
        print('\n'+try_catch(self.controller.get_filter_help, [line.strip()])+'\n')

## Output info to console

    def help_barfgenegff(self):
        print("Usage: barfgenegff <gene_id>\n")
        print("Prints gff entry for corresponding gene to console.\n")

    def do_barfgenegff(self, line):
        print(try_catch(self.controller.barf_gene_gff, [line]))

    def help_barfseq(self):
        print("Usage: barfseq <seq_id> <start_index> <end_index>\n")
        print("Prints (sub)sequence to console.\n")

    def do_barfseq(self, line):
        print(try_catch(self.controller.barf_seq, [line]))

    def help_barfcdsseq(self):
        print("Usage: barfcdsseq <mrna>\n")
        print("Prints CDS's whole sequence\n")

    def do_barfcdsseq(self, line):
        print(try_catch(self.controller.barf_cds_seq, [line]))

    def help_barfgenetbl(self):
        print("TODO")   # TODO

    def do_barfgenetbl(self, line):
        print(try_catch(self.controller.barf_gene_tbl, [line]))

    def help_info(self):
        print("Usage: info\n")
        print("Prints summary statistics about original genome (from file)" +\
                " and modified genome (filters applied). May take a moment to run.")

    def do_info(self, line):
        print(try_catch(self.controller.stats, None))

## Output info to file

    def help_writetbl(self):
        print("Usage: writetbl <file_name>\n")
        print("Write a sweet feature table to the specified file.\n")

    def do_writetbl(self, line):
        print(try_catch(self.controller.write_tbl, [line]))


##############################################

class LoadCmd(cmd.Cmd):

    intro = "Welcome to the GAG LOAD menu.\n"+\
            "Type the path to a folder containing your .fasta and .gff files.\n"+\
            "To use the current directory, just hit enter.\n"+\
            "You can type 'home' at any time to return to the main GAG console.\n"+\
            "You'll be returned automatically once your genome is loaded.\n"

    def __init__(self, prompt_prefix, controller, path_to_load):
        cmd.Cmd.__init__(self)
        self.prompt = prompt_prefix[:-2] + " LOAD> "
        self.controller = controller
        if path_to_load:
            self.cmdqueue = [path_to_load] # Execute default method with path as arg
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_home(self):
        print("Exit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True

    def exit_if_genome_loaded(self):
        if self.controller.seqs:
            return True

    def help_load(self):
        print(intro)

    def emptyline(self):
        self.default(".")
        return self.exit_if_genome_loaded()

    def default(self, line):
        if self.controller.seqs:
            print("Clearing genome ...")
            self.controller.clear_seqs()
        try_catch(self.controller.load_folder, [line])
        return self.exit_if_genome_loaded()

################################################

if __name__ == '__main__':
    GagCmd().cmdloop()
