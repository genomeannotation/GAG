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

    intro = "                                                                                \n                                        ,,                                      \n              ~++=:,                  =OMNZ:                  :=?77I=,          \n           :7DMMMMMMD7:              ~DMMMMZ               =$DMMMMMMMM8+        \n         :7NMMMMMMMMMMN7,            $MMMMMN=            ~ZMMMMMMMMMMMMZ,       \n        ?NMMMNZ+~+IONMMN?           INMMNMMMI           INMMMDZ?~,,:+I?:        \n       +NMMN7:      :?$I:          INMMNZMMMZ         :$MMMD?                   \n      :8MMM?                      7MMMD+=NMM8,       :OMMMO~                    \n      7MMMZ,                     7MMMD= ,8MMM?       ZMMMZ,                     \n     :8MMN+                    ,$MMMN+   IMMMO      ~NMMD~     ,+$ZI+==~,       \n     ~DMMD~    ~7OZ$Z8NNZ:    +DMMMD$I?~:~DMMN+     +NMM8      IMMMMMMMMO:      \n      $MMM8:  ,OMMMMMMMMN=   INMMN$$MMMMMMMMMMZ     ~DMMMI     ~8MMMMMMMZ:      \n      ,ZMMMD?, ?NMMMMMNZ~   INMMD+ =ONMMMMMMMMN~     +DMMM$:     ,=$MMMD:       \n       ,$MMMMM87+IONMMMI   ?NMMN=     :=+??OMMM$      =8MMMM87=, ,+8MMMZ        \n         +8MMMMMMMMMMMO~  ~DMMN?           =NMM$       :$NMMMMMMNMMMMMO:        \n           :I8NMMMMNO+,   :OMD?             +Z7~         :I8NMMMMMMMDI,         \n               :~~:        ,,,                               ~?7$7?~            \n                                                                                \n                                                                                "+\
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

    def help_load(self):
        print("This command takes you the GAG LOAD menu. There you can specify the location of")
        print("your files and load them into memory.")
        print("Alternately, just type 'load <path>' and avoid the submenu altogether.\n")

    def do_load(self, line):
        path_to_load = line.strip()
        loadcmd = LoadCmd(self.prompt, self.controller, path_to_load)
        loadcmd.cmdloop()

    def help_write(self):
        print("This command takes you to the GAG WRITE menu. There you can write genomic data")
        print("to the screen or to a file. You can write at the CDS, mRNA, gene, sequence or genome level.")
        print("Available formats: fasta, gff, tbl.\n")
        
    def do_write(self, line):
        # TODO confirm genome loaded
        writecmd = WriteCmd(self.prompt, self.controller, line) # Pass args to next console for parsing
        writecmd.cmdloop()

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

    helptext = "This is the GAG LOAD menu.\n"+\
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
        print("Exit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True

    def genome_loaded(self):
        if self.controller.seqs:
            return True

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

    helptext = "Welcome to the GAG WRITE menu.\n"+\
            "You can write at the cds, gene, seq or genome level. Please type your choice.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n"+\
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
        print("Exit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True

    def help_write(self):
        print(self.helptext)

    def emptyline(self):
        print(self.helptext)

    def do_fasta(self, line):
        fastacmd = WriteFastaCmd(self.prompt, self.controller, self.context, line)
        fastacmd.cmdloop()
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

class WriteFastaCmd(cmd.Cmd):

    helptext = "Welcome to the GAG WRITE CDS menu.\n"+\
            "You can write a CDS to fasta, gff or tbl format.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n"+\
            "fasta, gff or tbl?"

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
        print("Exit this console and return to the main GAG console.\n")

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

    helptext = "Welcome to the GAG WRITE CDS FASTA menu.\n"+\
            "Please type the mRNA id that corresponds to the CDS you want to write.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n"

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
        print("Exit this console and return to the main GAG console.\n")

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

    helptext = "Welcome to the GAG WRITE GENE menu.\n"+\
            "You can write a gene to gff or tbl format.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n"+\
            "gff or tbl?"

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
        print("Exit this console and return to the main GAG console.\n")

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

    helptext = "Welcome to the GAG WRITE GENE GFF menu.\n"+\
            "Please type the gene id that you want to write.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n"

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
        print("Exit this console and return to the main GAG console.\n")

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

    helptext = "Welcome to the GAG WRITE GENE TBL menu.\n"+\
            "Please type the gene id that you want to write.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n"

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
        print("Exit this console and return to the main GAG console.\n")

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

    helptext = "Welcome to the GAG WRITE SEQ menu.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n"+\
            "Please type the seq id you wish to write, followed by the start and stop bases.\n"

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
        print("Exit this console and return to the main GAG console.\n")

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

    helptext = "Welcome to the GAG WRITE GENOME menu.\n"+\
            "You can write a genome to fasta, gff or tbl file,\n"+\
            "or you can write all three files to a folder.\n"+\
            "(Type 'home' at any time to return to the main GAG console.)\n"+\
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
        print("Exit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def do_fasta(self, line):
        print("Genome to fasta coming soon!")

    def do_gff(self, line):
        print("Genome to gff coming soon!")

    def do_tbl(self, line):
        print("Genome to tbl coming soon!")

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
