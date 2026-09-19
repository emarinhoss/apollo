import os
import cmd
import readline

class WxInteractivePlot(cmd.Cmd):

    def __init__(self, options):
        self.options = options
        cmd.Cmd.__init__(self)
        self.prompt = "plot> "
        self.vars = {'inputFile' : None,
                     'frame' : None,
                     'variable' : None,
                     'component' : None,
                     'transformVariable' : None,
                     'transformsFile' : None}

    ## Command definitions ##
    def do_hist(self, args):
        """Print a list of commands that have been entered"""
        print(self._hist)
    def do_exit(self, args):
        """Exits from the console"""
        return -1

    def do_list(self, args):
        """Show list of variables which can be set"""
        print(self.vars.keys())

    def do_show(self, args):
        """Show value of option"""
        print(self.vars[args])

    def do_set(self, args):
        """Set value of variable"""
        lhs, rhs = args.split('=')
        lhs = lhs.strip()
        rhs = rhs.strip()
        self.vars[lhs] = rhs

    ## Command definitions to support Cmd object functionality ##
    def do_EOF(self, args):
        """Exit on system end of file character"""
        return self.do_exit(args)

    def do_shell(self, args):
        """Pass command to a system shell when line begins with '!'"""
        os.system(args)

    def do_help(self, args):
        """Get help on commands
           'help' or '?' with no arguments prints a list of commands for which help is available
           'help <command>' or '? <command>' gives help on <command>
        """
        ## The only reason to define this method is for the help text in the doc string
        cmd.Cmd.do_help(self, args)

    ## Override methods in Cmd object ##
    def preloop(self):
        """Initialization before prompting user for commands.
           Despite the claims in the Cmd documentaion, Cmd.preloop() is not a stub.
        """
        cmd.Cmd.preloop(self)   ## sets up command completion
        self._hist    = []      ## No history yet
        self._locals  = {}      ## Initialize execution namespace for user
        self._globals = {}

    def postloop(self):
        """Take care of any unfinished business.
           Despite the claims in the Cmd documentaion, Cmd.postloop() is not a stub.
        """
        cmd.Cmd.postloop(self)   ## Clean up command completion
        print("Exiting...")

    def precmd(self, line):
        """ This method is called after the line has been input but before
            it has been interpreted. If you want to modifdy the input line
            before execution (for example, variable substitution) do it here.
        """
        self._hist += [ line.strip() ]
        return line

    def postcmd(self, stop, line):
        """If you want to stop the console, return something that evaluates to true.
           If you want to do some post command processing, do it here.
        """
        return stop

    def emptyline(self):    
        """Do nothing on empty input line"""
        pass

    def default(self, line):
        """Called on an input line when the command prefix is not recognized.
           In that case we execute the line as Python code.
        """
        try:
            exec(line, self._locals, self._globals)
        except Exception as e:
            print(e.__class__, ":", e)

if __name__ == '__main__':
        console = WxInteractivePlot()
        console.cmdloop()
