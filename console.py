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

    intro = "                                                                                \n                                        ,,                                      \n              ~++=:,                  =OMNZ:                  :=?77I=,          \n           :7DMMMMMMD7:              ~DMMMMZ               =$DMMMMMMMM8+        \n         :7NMMMMMMMMMMN7,            $MMMMMN=            ~ZMMMMMMMMMMMMZ,       \n        ?NMMMNZ+~+IONMMN?           INMMNMMMI           INMMMDZ?~,,:+I?:        \n       +NMMN7:      :?$I:          INMMNZMMMZ         :$MMMD?                   \n      :8MMM?                      7MMMD+=NMM8,       :OMMMO~                    \n      7MMMZ,                     7MMMD= ,8MMM?       ZMMMZ,                     \n     :8MMN+                    ,$MMMN+   IMMMO      ~NMMD~     ,+$ZI+==~,       \n     ~DMMD~    ~7OZ$Z8NNZ:    +DMMMD$I?~:~DMMN+     +NMM8      IMMMMMMMMO:      \n      $MMM8:   ,OMMMMMMMN=   INMMN$$MMMMMMMMMMZ     ~DMMMI     ~8MMMMMMMZ:      \n      ,ZMMMD?,   ?NMMMNZ~   INMMD+ =ONMMMMMMMMN~     +DMMM$:     ,=$MMMD:       \n       ,$MMMMM87+IONMMMI   ?NMMN=     :=+??OMMM$      =8MMMM87=, ,+8MMMZ        \n         +8MMMMMMMMMMMO~  ~D0GE?           =NMM$       :$NMMMMMMNMMMMMO:        \n           :I8NMMMMNO+,   :OMD?             +Z7~         :I8NMMMMMMMDI,         \n               :~~:        ,,,                               ~?7$7?~            \n                                                                                \n                                                                                "+\
    "\nWelcome to the GAG console!\n"+\
    "Type 'help' for available commands.\n"

    no_genome_message = "\nIt looks like no genome is currently loaded. Try the 'load' command.\n"+\
            "Type 'help load' to learn how to use it, or just 'help' for general advice.\n"

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

    def help_load(self):
        print("\nThis command takes you the GAG LOAD menu. There you can specify the location of")
        print("your files and load them into memory.")
        print("Alternately, just type 'load <path>' and avoid the submenu altogether.\n")

    def do_load(self, line):
        if self.controller.genome_is_loaded():
            path_to_load = line.strip()
            loadcmd = LoadCmd(self.prompt, self.controller, path_to_load)
            loadcmd.cmdloop()
        else:
            print(self.no_genome_message)

    def help_fix(self):
        print("\nThis command takes you to the GAG FIX menu. There you can apply fixes to the genome,")
        print("to resolve issues such as internal stops, terminal Ns, etc.")
        print("Fixes are applied when you type the 'info' command or when you write")
        print("the genome to a file.")
        print("Alternately, just type 'fix <name_of_fix> if you've done this before :)\n")

    def do_fix(self, line):
        if self.controller.genome_is_loaded():
            name_of_fix = line.strip()
            fixcmd = FixCmd(self.prompt, self.controller, name_of_fix)
            fixcmd.cmdloop()
        else:
            print(self.no_genome_message)

    def help_write(self):
        print("\nThis command takes you to the GAG WRITE menu. There you can write genomic data")
        print("to the screen or to a file. You can write at the CDS, mRNA, gene, sequence or genome level.")
        print("Available formats: fasta, gff, tbl.\n")
        
    def do_write(self, line):
        if self.controller.genome_is_loaded():
            writecmd = WriteCmd(self.prompt, self.controller, line) # Pass args to next console for parsing
            writecmd.cmdloop()
        else:
            print(self.no_genome_message)

    def help_exit(self):
        print("\nExit this console.\n")

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

    def help_info(self):
        print("\nPrints summary statistics about original genome (from file)" +\
                " and modified genome (filters and fixes applied).")
        print("May take a moment to run.\n")

    def do_info(self, line):
        if self.controller.genome_is_loaded():
            print(try_catch(self.controller.stats, None))
        else:
            print(self.no_genome_message)


##############################################

class FixCmd(cmd.Cmd):

    helptext = "\nThis is the GAG FIX menu.\n"+\
            "You can apply the following fixes: "+\
            "terminal_ns, internal_stops, start_stop_codons.\n"+\
            "(You can type 'home' at any time to return to the main GAG console.)\n"+\
            "Type the name of a fix to enable it. Type 'all' to enable everything.\n\n"+\
            "terminal_ns, internal_stops, start_stop_codons or all?\n"

    def __init__(self, prompt_prefix, controller, name_of_fix):
        cmd.Cmd.__init__(self)
        self.prompt = prompt_prefix[:-2] + " FIX> "
        self.controller = controller
        if name_of_fix:
            self.cmdqueue = [name_of_fix] # Execute default method with path as arg
        else:
            print(self.helptext)
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True

    def help_terminal_ns(self):
        print("\nWith this fix enabled, GAG will inspect the beginning and end of each")
        print("sequence for the presence of unknown bases ('N' or 'n'). If found, they")
        print("will be removed. Indices of features on the sequence will be automatically")
        print("adjusted, and any mRNA that extended into the trimmed region will be removed.\n")

    def do_terminal_ns(self, line):
        self.controller.fix_terminal_ns()

    def help_internal_stops(self):
        print("\nEnabling this fix causes GAG to inspect each CDS for the presence of internal stops.")
        print("If any are found, a six-frame translation is performed. If the translation in another")
        print("phase is found to have no internal stops, this phase is considered to be the correct one")
        print("and the CDS is adjusted accordingly. If all translations contain internal stops, the CDS")
        print("and its parent mRNA are discarded.\n")

    def do_internal_stops(self, line):
        self.controller.fix_internal_stops()

    def help_start_stop_codons(self):
        print("\nSelecting this fix will cause GAG to inspect the first and last three bases of each CDS")
        print("to determine if it contains a valid start and/or stop codon. Information in the original")
        print("GFF regarding start_codon or stop_codon features is disregarded.\n")

    def do_start_stop_codons(self, line):
        self.controller.fix_start_stop_codons()

    def help_all(self):
        print("\nEnables all available fixes.\n")

    def do_all(self, line):
        self.controller.fix_terminal_ns()
        self.controller.fix_internal_stops()
        self.controller.fix_start_stop_codons()

    def emptyline(self):
        print(self.helptext)

    def default(self, line):
        print("Sorry, can't fix " + line)
        print(self.helptext)

##############################################

class LoadCmd(cmd.Cmd):

    helptext = "\nThis is the GAG LOAD menu.\n"+\
            "Type the path to a folder containing your .fasta and .gff files.\n"+\
            "To use the current directory, just hit enter.\n"+\
            "You can type 'home' at any time to return to the main GAG console.\n"+\
            "You'll be returned automatically once your genome is loaded.\n\n"+\
            "Folder path?\n"

    def __init__(self, prompt_prefix, controller, path_to_load):
        cmd.Cmd.__init__(self)
        self.prompt = prompt_prefix[:-2] + " LOAD> "
        self.controller = controller
        if path_to_load:
            self.cmdqueue = [path_to_load] # Execute default method with path as arg
        else:
            print(self.helptext)
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True

    def genome_loaded(self):
        return self.controller.genome_is_loaded()

    def help_load(self):
        print(self.helptext)

    def emptyline(self):
        self.default(".")
        return self.genome_loaded()

    def default(self, line):
        if self.controller.seqs:
            print("Clearing genome ...")
            self.controller.clear_seqs()
        try_catch(self.controller.load_folder, [line])
        if self.genome_loaded():
            return True
        else:
            print(self.helptext)

################################################

class WriteCmd(cmd.Cmd):

    helptext = "\nWelcome to the GAG WRITE menu.\n"+\
            "You can write at the cds, gene, seq or genome level. Please type your choice.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
            "cds, gene, seq or genome?\n"

    def __init__(self, prompt_prefix, controller, line):
        cmd.Cmd.__init__(self)
        self.prompt = prompt_prefix[:-2] + " WRITE> "
        self.controller = controller
        self.context = {"go_home": False}
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True

    def help_write(self):
        print(self.helptext)

    def emptyline(self):
        print(self.helptext)

    def do_cds(self, line):
        cdscmd = WriteCDSCmd(self.prompt, self.controller, self.context, line)
        cdscmd.cmdloop()
        if self.context["go_home"]:
            return True

    def do_gene(self, line):
        genecmd = WriteGeneCmd(self.prompt, self.controller, self.context, line)
        genecmd.cmdloop()
        if self.context["go_home"]:
            return True

    def do_seq(self, line):
        seqcmd = WriteSeqCmd(self.prompt, self.controller, self.context, line)
        seqcmd.cmdloop()
        if self.context["go_home"]:
            return True

    def do_genome(self, line):
        genomecmd = WriteGenomeCmd(self.prompt, self.controller, self.context, line)
        genomecmd.cmdloop()
        if self.context["go_home"]:
            return True

    def default(self, line):
        response = "Sorry, I don't know how to write " + line + "."
        response += "Please choose 'cds', 'gene', 'seq' or 'genome',"
        response += "or type 'home' to return to the main menu."
        print(response)

################################################

class WriteCDSCmd(cmd.Cmd):

    helptext = "\nWelcome to the GAG WRITE CDS menu.\n"+\
            "You can write a CDS to fasta, gff or tbl format.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
            "fasta, gff or tbl?\n"

    def __init__(self, prompt_prefix, controller, context, line):
        cmd.Cmd.__init__(self)
        self.prompt = prompt_prefix[:-2] + " CDS> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def do_fasta(self, line):
        cdsfastacmd = WriteCDSFastaCmd(self.prompt, self.controller, self.context, line)
        cdsfastacmd.cmdloop()
        if self.context["go_home"]:
            return True

    def do_gff(self, line):
        print("CDS to gff coming soon!")

    def do_tbl(self, line):
        print("CDs to tbl coming soon!")
    
    def help_writecds(self):
        print(self.helptext)

    def emptyline(self):
        print(self.helptext)

    def default(self, line):
        print("Sorry, I don't know how to write to " + line + " format.")
        print(self.helptext)

################################################

class WriteCDSFastaCmd(cmd.Cmd):

    helptext = "\nWelcome to the GAG WRITE CDS FASTA menu.\n"+\
            "Please type the mRNA id that corresponds to the CDS you want to write.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
            "mRNA id?\n"

    def __init__(self, prompt_prefix, controller, context, line):
        cmd.Cmd.__init__(self)
        self.prompt = prompt_prefix[:-2] + " FASTA> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def help_writecdsfasta(self):
        print(self.helptext)

    def emptyline(self):
        print(self.helptext)

    def default(self, line):
        print(try_catch(self.controller.barf_cds_seq, [line]))
        # TODO try mrna id; if not, give samples
        # TODO allow write to file
        # TODO return home if successful
        pass

################################################

class WriteGeneCmd(cmd.Cmd):

    helptext = "\nWelcome to the GAG WRITE GENE menu.\n"+\
            "You can write a gene to gff or tbl format.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
            "gff or tbl?\n"

    def __init__(self, prompt_prefix, controller, context, line):
        cmd.Cmd.__init__(self)
        self.prompt = prompt_prefix[:-2] + " GENE> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def do_gff(self, line):
        genegffcmd = WriteGeneGFFCmd(self.prompt, self.controller, self.context, line)
        genegffcmd.cmdloop()
        if self.context["go_home"]:
            return True

    def do_tbl(self, line):
        genetblcmd = WriteGeneTBLCmd(self.prompt, self.controller, self.context, line)
        genetblcmd.cmdloop()
        if self.context["go_home"]:
            return True

    def do_fasta(self, line):
        cdsfastacmd = WriteGeneFastaCmd(self.prompt, self.controller, self.context, line)
        cdsfastacmd.cmdloop()
        if self.context["go_home"]:
            return True

    def help_writegene(self):
        print(self.helptext)

    def emptyline(self):
        print(self.helptext)

    def default(self, line):
        print("Sorry, I don't know how to write to " + line + " format.")
        print(self.helptext)

################################################

class WriteGeneGFFCmd(cmd.Cmd):

    helptext = "\nWelcome to the GAG WRITE GENE GFF menu.\n"+\
            "Please type the gene id that you want to write.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
            "gene id?\n"

    def __init__(self, prompt_prefix, controller, context, line):
        cmd.Cmd.__init__(self)
        self.prompt = prompt_prefix[:-2] + " GFF> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def help_writegenegff(self):
        print(self.helptext)

    def emptyline(self):
        print(self.helptext)

    def default(self, line):
        print(try_catch(self.controller.barf_gene_gff, [line]))
        # TODO try gene id; if not, give samples
        # TODO allow write to file
        # TODO return home if successful
        pass

################################################

class WriteGeneTBLCmd(cmd.Cmd):

    helptext = "\nWelcome to the GAG WRITE GENE TBL menu.\n"+\
            "Please type the gene id that you want to write.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
            "gene id?\n"

    def __init__(self, prompt_prefix, controller, context, line):
        cmd.Cmd.__init__(self)
        self.prompt = prompt_prefix[:-2] + " TBL> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def help_writegenegff(self):
        print(self.helptext)

    def emptyline(self):
        print(self.helptext)

    def default(self, line):
        print(try_catch(self.controller.barf_gene_tbl, [line]))
        # TODO try gene id; if not, give samples
        # TODO allow write to file
        # TODO return home if successful
        pass

################################################

class WriteSeqCmd(cmd.Cmd):

    helptext = "\nWelcome to the GAG WRITE SEQ menu.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n"+\
            "Please type the seq id you wish to write, followed by the start and stop bases.\n\n"+\
            "seq id [start base] [stop base]?\n"

    def __init__(self, prompt_prefix, controller, context, line):
        cmd.Cmd.__init__(self)
        self.prompt = prompt_prefix[:-2] + " SEQ> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def help_writeseq(self):
        print(self.helptext)

    def emptyline(self):
        print(self.helptext)

    def default(self, line):
        print(try_catch(self.controller.barf_seq, [line]))
        # TODO make start/stop base optional in consolecontroller method
        # TODO try seq; if not provide hints
        # TODO If failed print help message

################################################

class WriteGenomeCmd(cmd.Cmd):

    helptext = "\nWelcome to the GAG WRITE GENOME menu.\n"+\
            "You can write a genome to fasta, gff or tbl file,\n"+\
            "or you can write all three files to a folder.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
            "fasta, gff, tbl, or all?"

    def __init__(self, prompt_prefix, controller, context, line):
        cmd.Cmd.__init__(self)
        self.prompt = prompt_prefix[:-2] + " GENOME> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return cmd.Cmd.precmd(self, line)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def do_fasta(self, line):
        print("Genome to fasta coming soon!")

    def do_gff(self, line):
        print("Genome to gff coming soon!")

    def do_tbl(self, line):
        # TODO verify that line is valid path first
        # TODO return home if successful
        print(try_catch(self.controller.write_tbl, [line]))

    def do_all(self, line):
        # TODO verify line is valid path first? or does console controller do that?
        # TODO return home if successful
        print(try_catch(self.controller.barf_folder, [line]))
    
    def help_writecds(self):
        print(self.helptext)

    def emptyline(self):
        print(self.helptext)

    def default(self, line):
        print("Sorry, I don't know how to write to " + line + " format.")
        print(self.helptext)

################################################

if __name__ == '__main__':
    GagCmd().cmdloop()
