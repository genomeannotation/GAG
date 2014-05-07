#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import ast
import cmd
import readline
import sys
import traceback
import types
from src.console_controller import ConsoleController

def try_catch(command, args=None):
    try:
        if args is not None:
            return command(*args)
        else:
            return command()
    except:
        print("Sorry, that command raised an exception. Here's what I know:\n")
        print(traceback.format_exc())
        
###################################################################################################
# Define a custom cmd line base class to fix some bugs in cmd.Cmd
###################################################################################################

# TODO This class introduces some cool opportunities - we might be able to put all the duplicated
# functionality of subconsoles (i.e. home related code) into this class so the subconsoles' code
# can be more concise

class GagCmdBase(cmd.Cmd):

    def __init__(self):
        cmd.Cmd.__init__(self)
        self.helptext = ""
        readline.set_history_length(1000)
        try:
            readline.read_history_file('.gaghistory')
        except IOError:
            sys.stderr.write("No .gaghistory file available...\n")
        
    def get_names(self):
        return dir(self)

    def precmd(self, line):
        readline.write_history_file('.gaghistory')
        return line

    def emptyline(self):
        print(self.helptext)

###################################################################################################
# End cmd base class
###################################################################################################

class GagCmd(GagCmdBase):

## Setup, loading and saving sessions, exit, etc.

    intro = "                                                                                \n                                        ,,                                      \n              ~++=:,                  =OMNZ:                  :=?77I=,          \n           :7DMMMMMMD7:              ~DMMMMZ               =$DMMMMMMMM8+        \n         :7NMMMMMMMMMMN7,            $MMMMMN=            ~ZMMMMMMMMMMMMZ,       \n        ?NMMMNZ+~+IONMMN?           INMMNMMMI           INMMMDZ?~,,:+I?:        \n       +NMMN7:      :?$I:          INMMNZMMMZ         :$MMMD?                   \n      :8MMM?                      7MMMD+=NMM8,       :OMMMO~                    \n      7MMMZ,                     7MMMD= ,8MMM?       ZMMMZ,                     \n     :8MMN+                    ,$MMMN+   IMMMO      ~NMMD~     ,+$ZI+==~,       \n     ~DMMD~    ~7OZ$Z8NNZ:    +DMMMD$I?~:~DMMN+     +NMM8      IMMMMMMMMO:      \n      $MMM8:   ,OMMMMMMMN=   INMMN$$MMMMMMMMMMZ     ~DMMMI     ~8MMMMMMMZ:      \n      ,ZMMMD?,   ?NMMMNZ~   INMMD+ =ONMMMMMMMMN~     +DMMM$:     ,=$MMMD:       \n       ,$MMMMM87+IONMMMI   ?NMMN=     :=+??OMMM$      =8MMMM87=, ,+8MMMZ        \n         +8MMMMMMMMMMMO~  ~D0GE?           =NMM$       :$NMMMMMMNMMMMMO:        \n           :I8NMMMMNO+,   :OMD?             +Z7~         :I8NMMMMMMMDI,         \n               :~~:        ,,,                               ~?7$7?~            \n                                                                                \n                                                                                "+\
    "\nWelcome to the GAG console!\n"+\
    "Type 'help' for available commands.\n"

    no_genome_message = "\nIt looks like no genome is currently loaded. Try the 'load' command.\n"+\
            "Type 'help load' to learn how to use it, or just 'help' for general advice.\n"

    def __init__(self):
        GagCmdBase.__init__(self)
        self.prompt = "GAG> "
        self.controller = ConsoleController() 

    def precmd(self, line):
        line = GagCmdBase.precmd(self, line)
        print('----------------------------------------')
        return line
        
    def postcmd(self, stop, line):
        stop = GagCmdBase.postcmd(self, stop, line)
        print('----------------------------------------')
        return stop

    def help_load(self):
        print("\nThis command takes you the GAG LOAD menu. There you can specify the location of")
        print("your files and load them into memory.")
        print("Alternately, just type 'load <path>' and avoid the submenu altogether.\n")

    def do_load(self, line):
        path_to_load = line.strip()
        loadcmd = LoadCmd(self.prompt, self.controller, path_to_load)
        loadcmd.cmdloop()

    def help_annotate(self):
        print("\nThis command takes you to the GAG ANNOTATE menu. There you can specify the location of")
        print("a table of annotations. GAG will add these annotations to your genome.")
        print("Alternately, just type 'annotate <path_to_table> and avoid the submenu altogether.\n")

    def do_annotate(self, line):
        if self.controller.genome_is_loaded():
            annocmd = AnnoCmd(self.prompt, self.controller, line)
            annocmd.cmdloop()
        else:
            print(self.no_genome_message)

    def help_flag(self):
        print("\nThis command takes you to the GAG FLAG menu. There you can flag features based")
        print("on certain criteria to mark suspicious data. Alternately, just type:\n")
        print("'flag <criteria> <value>'\n")
        print("if you've done this before :)\n")

    def do_flag(self, line):
        if self.controller.genome_is_loaded():
            filtercmd = FlagCmd(self.prompt, self.controller, line)
            filtercmd.cmdloop()
        else:
            print(self.no_genome_message)

    def help_trim(self):
        print("\nThis command takes you to the GAG TRIM menu, where you can ")
        print("trim regions of the genome using a .bed file.\n")
        print("Alternately, just type 'trim <path_to_file.bed>'\n")

    def do_trim(self, line):
        if self.controller.genome_is_loaded():
            trimcmd = TrimCmd(self.prompt, self.controller, line)
            trimcmd.cmdloop()
        else:
            print(self.no_genome_message)
    
    def help_remove(self):
        print("\nThis command takes you to the GAG REMOVE menu. There you can remove features based")
        print("on certain criteria to filter out bad data. Alternately, just type:\n")
        print("'remove <criteria> <value>'\n")
        print("if you've done this before :)\n")

    def do_remove(self, line):
        if self.controller.genome_is_loaded():
            filtercmd = RemoveCmd(self.prompt, self.controller, line)
            filtercmd.cmdloop()
        else:
            print(self.no_genome_message)
            
    def help_fix(self):
        print("\nThis command takes you to the GAG FIX menu. There you can apply fixes to the genome,")
        print("to resolve issues such as terminal Ns, invalid start and stop codons, etc.")
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
        print("\nThis command takes you to the GAG WRITE menu.")
        print("There you can write your genome to a folder.\n")
        
    def do_write(self, line):
        if self.controller.genome_is_loaded():
            writecmd = WriteCmd(self.prompt, self.controller, line) # Pass args to next console for parsing
            writecmd.cmdloop()
        else:
            print(self.no_genome_message)

    def help_view(self):
        print("\nThis command takes you to the GAG VIEW menu.")
        print("There you can view cds, gene or sequence features on screen.\n")

    def do_view(self, line):
        if self.controller.genome_is_loaded():
            viewcmd = ViewCmd(self.prompt, self.controller, line)
            viewcmd.cmdloop()
        else:
            print(self.no_genome_message)

    def help_exit(self):
        print("\nExit this console.\n")

    def do_exit(self, line):
        return True
        
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

class FlagCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nThis is the GAG FLAG menu.\n"+\
        "(You can type 'home' at any time to return to the main GAG console.)\n\n"+\
        " - Here you can flag features based on certain criteria.\n"+\
        " - When you flag features, they remain in the genome but are written with a 'gag_flagged'\n"+\
        "   annotation so they can be reviewed manually, in a genome browser for instance.\n"+\
        " - Typing info in the main GAG console will tell you the total amount of features that\n"+\
        "   have been flagged.\n"+\
        " - You can get a summary of the current flag criteria by typing 'summary'\n\n"+\
        "You can flag:\n\n"+\
        "\n".join(controller.filter_mgr.filters.keys())+ "\n"
        self.prompt = prompt_prefix[:-2] + " FLAG> "
        self.controller = controller
        self.context = {"go_home": False}
        if line:
            self.cmdqueue = [line] # Execute default method with path as arg
        else:
            print(self.helptext)
            
        # Set up filter arg do functions
        for filt_name in controller.filter_mgr.filters.keys():
            # First real closure #teddy'sallgrownup
            # traps arg variable
            def do_arg(slf, line, filter_name = filt_name):
                filtercmd = FilterArgCmd(slf.prompt, slf.controller, slf.context, line, 'FLAG', filter_name)
                filtercmd.cmdloop()
                if slf.context["go_home"]:
                    return True
            setattr(self, 'do_'+filt_name, types.MethodType(do_arg, self))

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True
    
    def do_summary(self, line):
        for filt_name, filt in self.controller.filter_mgr.filters.iteritems():
            if not filt.remove:
                print('flag ' + filt_name + ' ' + str(filt.arg))

    def default(self, line):
        print("Sorry, can't flag " + line)
        print(self.helptext)
        
##############################################

class TrimCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nThis is the GAG TRIM menu.\n"+\
        "(You can type 'home' at any time to return to the main GAG console.)\n"+\
        "Please type the path to a .bed file containing regions to trim.\n\n"
        self.prompt = prompt_prefix[:-2] + " TRIM> "
        self.controller = controller
        self.context = {"go_home": False}
        if line:
            self.cmdqueue = [line] # Execute default method with path as arg
        else:
            print(self.helptext)
            
    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True
    
    def default(self, line):
        try_catch(self.controller.trim_from_file, [line])
        return True
        
##############################################

class AnnoCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nThis is the GAG ANNOTATE menu.\n"+\
        "(You can type 'home' at any time to return to the main GAG console.)\n"+\
        "Please type the path to a table containing your annotations.\n\n"
        self.prompt = prompt_prefix[:-2] + " ANNOTATE> "
        self.controller = controller
        self.context = {"go_home": False}
        if line:
            self.cmdqueue = [line] # Execute default method with path as arg
        else:
            print(self.helptext)
            
    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True
    
    def default(self, line):
        try_catch(self.controller.annotate_from_file, [line])
        return True

##############################################
        
class RemoveCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nThis is the GAG REMOVE menu.\n"+\
        "(You can type 'home' at any time to return to the main GAG console.)\n\n"+\
        " - Here you can remove features based on certain criteria.\n"+\
        " - You can view the effects of removal by typing 'info' in the main GAG console.\n"+\
        " - You can get a summary of the current remove criteria by typing 'summary'\n\n"+\
        "You can flag:\n\n"+\
        "You can remove:\n\n"+\
        "\n".join(controller.filter_mgr.filters.keys())+ "\n"+\
        "\nYou can also remove sequences, genes and mRNAs from a file by typing 'from_file'.\n\n"
        self.prompt = prompt_prefix[:-2] + " REMOVE> "
        self.controller = controller
        self.context = {"go_home": False}
        if line:
            self.cmdqueue = [line] # Execute default method with path as arg
        else:
            print(self.helptext)
            
        # Set up filter arg do functions
        for filt_name in controller.filter_mgr.filters.keys():
            # First real closure #teddy'sallgrownup
            # traps arg variable
            def do_arg(slf, line, filter_name = filt_name):
                filtercmd = FilterArgCmd(slf.prompt, slf.controller, slf.context, line, 'REMOVE', filter_name)
                filtercmd.cmdloop()
                if slf.context["go_home"]:
                    return True
            setattr(self, 'do_'+filt_name, types.MethodType(do_arg, self))

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True
    
    def do_summary(self, line):
        for filt_name, filt in self.controller.filter_mgr.filters.iteritems():
            if filt.remove and filt.arg != 0:
                print('remove ' + filt_name + ' ' + str(filt.arg))

    def do_from_file(self, line):
        if line:
            print("")
            try_catch(self.controller.remove_from_file, [line])
            print("")
        else:
            line = raw_input("Path to file? ")
            print("")
            try_catch(self.controller.remove_from_file, [line])
            print("")
        return True

    def default(self, line):
        print("Sorry, can't remove " + line)
        print(self.helptext)
        
##############################################

# Generic filter arg setter/getter command console
class FilterArgCmd(GagCmdBase):

    # Note: filter_mode is either 'FLAG' or 'REMOVE'
    def __init__(self, prompt_prefix, controller, context, line, filter_mode, filter_name):
        GagCmdBase.__init__(self)
        self.helptext = "This is the "+filter_mode+" "+filter_name+" console. Type a value\n" \
                        "to set the criteria, or simply press enter to see the current value.\n"
        
        self.prompt = prompt_prefix[:-2] + ' ' + filter_name +'> '
        self.controller = controller
        self.context = context
        self.filter_mode = filter_mode
        self.filter_name = filter_name
        if line:
            self.cmdqueue = [line] # Execute default method with path as arg
        else:
            print(self.helptext)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context['go_home'] = True
        return True

    def emptyline(self):
        print(self.filter_name+": "+str(try_catch(self.controller.get_filter_arg, [self.filter_name]))+"\n\n")

    def default(self, line):
        line = line.strip()
        if line:
            # Validate the value
            got_type = 'invalid'
            try:
                got_type = type(ast.literal_eval(line)).__name__
            except:
                pass
            expected_type = type(self.controller.get_filter_arg(self.filter_name)).__name__
            # it's a number
            if got_type != expected_type:
                print("Failed to set "+self.filter_name+". Expected "+expected_type+".\n")
                if expected_type == 'bool':
                    print("A bool can be True or False (the caps matters).\n")
                elif expected_type == 'int':
                    if got_type == 'float': # Uber smart help
                        print("An int is any whole number. Try this: "+line.split('.')[0]+"\n")
                    else:
                        print("An int is any whole number (42 is valid, and 43.123 is invalid)\n")
                elif expected_type == 'float':
                    if got_type == 'int': # Uber smart help
                        print("A float is any decimal number. Try this: "+line+".0\n")
                    else:
                        print("A float is any decimal number (42.0, 43.123, -67.12 are all valid. If you want a whole number, it needs a .0 after it)\n")
                elif expected_type == 'str':
                    print("A str is any text. It is formatted \"like this\" or 'like this'\n")
                return False
            
            remove = False
            if self.filter_mode == 'REMOVE':
                remove = True

            try_catch(self.controller.apply_filter, [self.filter_name, line, remove])

            if self.filter_mode == 'REMOVE':
                print("\n"+self.filter_name+" "+line+" removed.\n")
            else: #TODO maybe throw an error if filter_mode isn't FLAG
                print("\n"+self.filter_name+" "+line+" flagged.\n")
            
            self.context['go_home'] = True
            return True

##############################################

class FixCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, name_of_fix):
        GagCmdBase.__init__(self)
        self.helptext = "\nThis is the GAG FIX menu.\n"+\
                "You can apply the following fixes: "+\
                "terminal_ns, start_stop_codons.\n"+\
                "(You can type 'home' at any time to return to the main GAG console.)\n"+\
                "Type the name of a fix to enable it. Type 'all' to enable everything.\n\n"+\
                "terminal_ns, start_stop_codons or all?\n"

        self.prompt = prompt_prefix[:-2] + " FIX> "
        self.controller = controller
        if name_of_fix:
            self.cmdqueue = [name_of_fix] # Execute default method with path as arg
        else:
            print(self.helptext)

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
        print("\n" + self.controller.fix_terminal_ns() + "\n")
        return True

    def help_start_stop_codons(self):
        print("\nSelecting this fix will cause GAG to inspect the first and last three bases of each CDS")
        print("to determine if it contains a valid start and/or stop codon. Information in the original")
        print("GFF regarding start_codon or stop_codon features is disregarded.\n")

    def do_start_stop_codons(self, line):
        print("\n" + self.controller.fix_start_stop_codons() + "\n")
        return True

    def help_all(self):
        print("\nEnables all available fixes.\n")

    def do_all(self, line):
        print("\n" + self.controller.fix_terminal_ns())
        print(self.controller.fix_start_stop_codons() + "\n")
        return True

    def default(self, line):
        print("Sorry, can't fix " + line)
        print(self.helptext)

##############################################

class LoadCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, path_to_load):
        GagCmdBase.__init__(self)
        self.helptext = "\nThis is the GAG LOAD menu.\n"+\
                "Type the path to a folder containing your .fasta and .gff files.\n"+\
                "To use the current directory, just hit enter.\n"+\
                "You can type 'home' at any time to return to the main GAG console.\n"+\
                "You'll be returned automatically once your genome is loaded.\n\n"+\
                "Folder path?\n"
        self.prompt = prompt_prefix[:-2] + " LOAD> "
        self.controller = controller
        if path_to_load:
            self.cmdqueue = [path_to_load] # Execute default method with path as arg
        else:
            print(self.helptext)

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
            print("\nGenome loaded. (Type 'info' to see summary statistics.)\n")
            return True
        else:
            print("\nError reading genome!")
            print(self.helptext)

################################################

class ViewCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nWelcome to the GAG VIEW menu.\n"+\
                "You can view at the cds, gene, or seq level. Please type your choice.\n"+\
                "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
                "cds, gene or seq?\n"
        self.prompt = prompt_prefix[:-2] + " VIEW> "
        self.controller = controller
        self.context = {"go_home": False}
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True

    def help_write(self):
        print(self.helptext)

    def do_cds(self, line):
        cdscmd = ViewCDSCmd(self.prompt, self.controller, self.context, line)
        cdscmd.cmdloop()
        if self.context["go_home"]:
            return True

    def do_gene(self, line):
        genecmd = ViewGeneCmd(self.prompt, self.controller, self.context, line)
        genecmd.cmdloop()
        if self.context["go_home"]:
            return True

    def do_seq(self, line):
        seqcmd = ViewSeqCmd(self.prompt, self.controller, self.context, line)
        seqcmd.cmdloop()
        if self.context["go_home"]:
            return True

    def do_genome(self, line):
        genomecmd = WriteGenomeCmd(self.prompt, self.controller, self.context, line)
        genomecmd.cmdloop()
        if self.context["go_home"]:
            return True

    def default(self, line):
        response = "\nSorry, I don't know how to display " + line + ".\n"
        response += "Please choose 'cds', 'gene' or 'seq',"
        response += "or type 'home' to return to the main menu.\n"
        print(response)
        print(self.helptext)

################################################

class ViewCDSCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, context, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nWelcome to the GAG VIEW CDS menu.\n"+\
                "You can view a CDS in fasta, gff or tbl format.\n"+\
                "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
                "fasta, gff or tbl?\n"
        self.prompt = prompt_prefix[:-2] + " CDS> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def do_fasta(self, line):
        cdsfastacmd = ViewCDSFastaCmd(self.prompt, self.controller, self.context, line)
        cdsfastacmd.cmdloop()
        if self.context["go_home"]:
            return True

    def do_gff(self, line):
        cdsgffcmd = ViewCDSGFFCmd(self.prompt, self.controller, self.context, line)
        cdsgffcmd.cmdloop()
        if self.context["go_home"]:
            return True

    def do_tbl(self, line):
        cdstblcmd = ViewCDSTBLCmd(self.prompt, self.controller, self.context, line)
        cdstblcmd.cmdloop()
        if self.context["go_home"]:
            return True
    
    def help_writecds(self):
        print(self.helptext)

    def default(self, line):
        print("Sorry, I don't know how to display in" + line + " format.")
        print(self.helptext)

################################################

class ViewCDSFastaCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, context, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nWelcome to the GAG VIEW CDS FASTA menu.\n"+\
                "Please type the mRNA id that corresponds to the CDS you want to view.\n"+\
                "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
                "mRNA id?\n"
        self.prompt = prompt_prefix[:-2] + " FASTA> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def help_viewcdsfasta(self):
        print(self.helptext)
        print(self.controller.get_n_mrna_ids(5))

    def emptyline(self):
        print(self.helptext)
        print(self.controller.get_n_mrna_ids(5))

    def default(self, line):
        mrna_id = line.strip()
        if self.controller.contains_mrna(mrna_id):
            print("\n" + try_catch(self.controller.barf_cds_seq, [line]) + "\n")
            self.context["go_home"] = True
            return True
        else:
            print("\nSorry, couldn't find mRNA id '" + mrna_id + "'.")
            print(self.controller.get_n_mrna_ids(5))

################################################

class ViewCDSGFFCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, context, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nWelcome to the GAG VIEW CDS GFF menu.\n"+\
                "Please type the mRNA id that corresponds to the CDS you want to view.\n"+\
                "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
                "mRNA id?\n"

        self.prompt = prompt_prefix[:-2] + " GFF> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def help_viewcdsgff(self):
        print(self.helptext)
        print(self.controller.get_n_mrna_ids(5))

    def emptyline(self):
        print(self.helptext)
        print(self.controller.get_n_mrna_ids(5))

    def default(self, line):
        mrna_id = line.strip()
        if self.controller.contains_mrna(mrna_id):
            print("\n" + try_catch(self.controller.cds_to_gff, [line]) + "\n")
            self.context["go_home"] = True
            return True
        else:
            print("\nSorry, couldn't find mRNA id '" + mrna_id + "'.")
            print(self.controller.get_n_mrna_ids(5))

################################################

class ViewCDSTBLCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, context, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nWelcome to the GAG VIEW CDS TBL menu.\n"+\
                "Please type the mRNA id that corresponds to the CDS you want to view.\n"+\
                "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
                "mRNA id?\n"
        self.prompt = prompt_prefix[:-2] + " TBL> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def help_viewcdstbl(self):
        print(self.helptext)
        print(self.controller.get_n_mrna_ids(5))

    def emptyline(self):
        print(self.helptext)
        print(self.controller.get_n_mrna_ids(5))

    def default(self, line):
        mrna_id = line.strip()
        if self.controller.contains_mrna(mrna_id):
            print("\n" + try_catch(self.controller.cds_to_tbl, [line]) + "\n")
            self.context["go_home"] = True
            return True
        else:
            print("\nSorry, couldn't find mRNA id '" + mrna_id + "'.")
            print(self.controller.get_n_mrna_ids(5))

################################################

class ViewGeneCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, context, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nWelcome to the GAG VIEW GENE menu.\n"+\
                "You can view a gene in gff or tbl format.\n"+\
                "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
                "gff or tbl?\n"
        self.prompt = prompt_prefix[:-2] + " GENE> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)

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

    def help_viewgene(self):
        print(self.helptext)

    def default(self, line):
        print("Sorry, I don't know how to display in " + line + " format.")
        print(self.helptext)

################################################

class WriteGeneGFFCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, context, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nWelcome to the GAG VIEW GENE GFF menu.\n"+\
                "Please type the gene id that you want to view.\n"+\
                "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
                "gene id?\n"
        self.prompt = prompt_prefix[:-2] + " GFF> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def help_viewgenegff(self):
        print(self.helptext)
        print(self.controller.get_n_gene_ids(5))

    def emptyline(self):
        print(self.helptext)
        print(self.controller.get_n_gene_ids(5))

    def default(self, line):
        gene_id = line.strip()
        if self.controller.contains_gene(gene_id):
            print("\n" + try_catch(self.controller.barf_gene_gff, [line]) + "\n")
            self.context["go_home"] = True
            return True
        else:
            print("\nSorry, couldn't find gene id '" + gene_id + "'.")
            print(self.controller.get_n_gene_ids(5))

################################################

class WriteGeneTBLCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, context, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nWelcome to the GAG VIEW GENE TBL menu.\n"+\
                "Please type the gene id that you want to view.\n"+\
                "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
                "gene id?\n"
        self.prompt = prompt_prefix[:-2] + " TBL> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def help_viewgenegff(self):
        print(self.helptext)
        print(self.controller.get_n_gene_ids(5))

    def emptyline(self):
        print(self.helptext)
        print(self.controller.get_n_gene_ids(5))

    def default(self, line):
        gene_id = line.strip()
        if self.controller.contains_gene(gene_id):
            print("\n" + try_catch(self.controller.barf_gene_tbl, [line]) + "\n")
            self.context["go_home"] = True
            return True
        else:
            print("\nSorry, couldn't find gene id '" + gene_id + "'.")
            print(self.controller.get_n_gene_ids(5))

################################################

class ViewSeqCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, context, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nWelcome to the GAG VIEW SEQ menu.\n"+\
                "(Type 'home' at any time to return to the main GAG console.)\n"+\
                "Please type the seq id you wish to view. To view a subsequence,\n"+\
                "follow the seq id with the start and stop bases.\n\n"+\
                "seq id [start base] [stop base]?\n"
        self.prompt = prompt_prefix[:-2] + " SEQ> "
        self.controller = controller
        self.context = context
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        self.context["go_home"] = True
        return True
    
    def help_viewseq(self):
        print(self.helptext)
        print(self.controller.get_n_seq_ids(5))

    def emptyline(self):
        print(self.helptext)
        print(self.controller.get_n_seq_ids(5))

    def default(self, line):
        args = line.split()
        seq_id = args[0]
        if self.controller.contains_seq(seq_id):
            print("\n" + try_catch(self.controller.barf_seq, [line]))
            self.context["go_home"] = True
            return True
        else:
            print("\nSorry, couldn't find seq id '" + seq_id + "'.")
            print(self.controller.get_n_seq_ids(5))
            print("seq id [start base] [stop base]?\n")

################################################

class WriteCmd(GagCmdBase):

    def __init__(self, prompt_prefix, controller, line):
        GagCmdBase.__init__(self)
        self.helptext = "\nWelcome to the GAG WRITE menu. From here you can write your genome to a folder.\n"+\
                "Please type a name for the folder to contain the files.\n"+\
                "The folder will be created -- in other words, don't give an existing folder.\n"+\
                "(Type 'home' at any time to return to the main GAG console.)\n\n"+\
                "folder name?\n"
        self.prompt = prompt_prefix[:-2] + " WRITE> "
        self.controller = controller
        if line:
            self.cmdqueue = [line] # Execute default method with passed-in line
        else:
            print(self.helptext)

    def help_home(self):
        print("\nExit this console and return to the main GAG console.\n")

    def do_home(self, line):
        return True
    
    def help_write(self):
        print(self.helptext)

    def default(self, line):
        if self.controller.can_write_to_path(line):
            print(try_catch(self.controller.barf_folder, [line]))
            return True
        else:
            print("\nSorry, can't write to " + line + "\n")
            print("folder name?\n")

################################################

if __name__ == '__main__':
    GagCmd().cmdloop()
