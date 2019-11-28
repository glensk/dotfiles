    #!/usr/bin/python


"""
my.py contains exclusive ONLY definitions of functions to use in other scripts. Just
>>> import my
and run a command with
>>> my.whateverdefintion
"""


def dos(infile=None):
    """
    create density of states
    """
    if infile is None:
        exit("Please define an inputfile")
    input = readfile(infile)
    print input


def run2(command=None):
    """
    constantly prints output, not just at the end
    """
    import subprocess
    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ''

    # Poll process for new output until finished
    for line in iter(process.stdout.readline, ""):
        print line,
        output += line
    process.wait()
    exitCode = process.returncode

    if (exitCode == 0):
        return output
    else:
        raise Exception(command, exitCode, output)


def run(befehl):
    """
    z.B. run('ls -l')
    z.B. run('frompath_stoich.sh')
    """
    import subprocess
    return subprocess.check_output(befehl, shell=True, stderr=subprocess.STDOUT)


def grep(string, list, options=None):
    if is_string(string) is not True:
        exit(str(string) + "is not a string")
    #    print "str:",string
    #    print "list:",list
    #    print "list",lis
    #    print "ll",len(list)
    #    print "ty",type(list)
    #    print "list?",is_list(list)
    if is_list(list) is not True:
        list = [list]
    #    print "list->",list
    import re

    if options == '-i':
        expr = re.compile(string, re.IGNORECASE)
    #        print "option -i"
    #        print ""
    else:
        expr = re.compile(string)
    out = []
    for text in list:
    #        print "jo11: ",text,"|||"+string+"|||"+"-->",re.search(string,text,re.IGNORECASE)
    #        print re.search(string,text)
    #        match = expr.search(text,re.IGNORECASE)
    #        print "match",match,";  expr:\""+string+"\";   text:\""+text+"\""
    #        print type(string),type(text)
        if options == '-i':
            match = re.search(string, text, re.IGNORECASE)
        else:
            match = re.search(string, text)
    #        print "MMM",match,text
        if match is not None:
            if options is None:
                out.insert(1, match.string)
                #print match.string
    #                print "hier0 search:",strint,match.string
            elif options == '-v':
                add = re.search(string, match.string).group()
                out.insert(1, add)
    #                print "hier1 search:",string,match.string
            elif options == '-i':
    #                print "hier2 search:",string,match.string
                add = text
                out.insert(1, add)
        else:
            pass
            #print None
    return out

    #    print "word:",word
    #    print "stringlist:",stringlist


    #    import re
    #    out=[]
    #    if len(stringlist) == 1:
    #
    #    for string in stringlist:
    #        print string
    #        if re.search(word, string) != None:
    #            add = re.search(word, string).group()
    #            out.insert(1,add)
    #    if len(out) != 0:
    #        return out[-1].replace('eV','')
    #    else:
    #        print "couldnt find number of kpoints from_path, I need exactlt one: I got",from_path_details_string().split("_")[:]
    #        my.exit()
    #    return out


def path_to_script(__file__):
    import os
    return os.path.realpath(__file__)

    #def get_definitions_of_script(inn):
    #    ## in sollte locals() sein
    #    definitions = []
    ##    in_ = str(locals()).split(",")
    #    in_ = str(inn).split(",")
    ##    print "IN:",in_
    #    for i in in_:
    #        ausdr = i.split("': <")
    ##        print "-------------->",ausdr,len(ausdr)
    #        if len(ausdr) == 2:
    #            if ausdr[1][0:8] == "function":
    #                out = ausdr[0].replace("{","").replace(" ","").replace("'","")
    #                definitions.append(out)
    #    return definitions

    #def get_definitions_of_script_bycatpath(path):
    #    out = []
    #    for line in readfile(path):
    #        if line != "": ## nur lines interessant welche nicht leer sind
    #        #            print "-----------> ",line.split()
    #        #            if line == "import argparse":
    #        #                my.exit()
    #            if line.split()[0] == "def":
    #                string = line.split()[1].split("(")[0]
    ##                print string
    #                out.append(string)
    #    return out


def get_definitions_of_current_script(script=None):
    #    print "XXXXXXXXXXX:script:",script
    if script is None:
    #        print ""
    #        print "--> ",which_script_show_allinspect()
    #        print ""
        path = which_scriptpath_isit_whew_i_was_called_from()
    else:
        path = checkfile(script)
    out = []
    #    print "PPPPP",path
    for line in readfile(path):
        if line is not "":  # nur lines interessant welche nicht leer sind
        #            print "-----------> ",line.split()
        #            if line == "import argparse":
        #                my.exit()
    #            print len(line),"line:::",line,":::______________"
    #            print "line.split",len(line.split())
            if len(line.split()) != 0:
                if line.split()[0] == "def":
                    string = line.split()[1].split("(")[0]
    #                print string
                    out.append(string)
    #    print "OUT:",out
    return out


def get_definitions_of_current_script_additionalinfo():
    """
    damit die defs funktionieren muss das entsprechende skript geladen werden: import test  oder impoer getdatabase.py
    docstrings muss sein: docstring = eval(deff).__doc__
    arguments muss sein: arguments = str(inspect.getcallargs(eval(deff)))
    """
    #    print "ok"
    path = which_scriptpath_isit_whew_i_was_called_from()
    #    print "KKK",path
    defs = get_definitions_of_current_script(script=path)

    #    #return path
    #
    #def defs_for_parser(defs):

    import inspect
    #callingscript = which_script_calledme()
    callingscript = which_scriptname_isit_whew_i_was_called_from()
    exec("import " + callingscript)  # statt import callingscript
    out = []
    #print "defs:",defs
    for deff in defs:
        """
        muss auch dafuer sorgen dass sich keine shortform widerholt!
        """
    #        print ""
    #        print "deff:",deff
    #        print "callstringscript:",callingscript
        docstring = eval(callingscript + "." + deff).__doc__
    #        print " "
    #        print "DEFF:",deff
    #        print "DOC>",docstring,"<DOC"

        ############### BEGIN DOCSTRING SHORTFORM ############
        ############### BEGIN DOCSTRING SHORTFORM ############
        if docstring is None:
            shortform = ""
            docstring = ""
        elif len(docstring) == 1:
            shortform = docstring
            docstring = ""
        else:
            shortform = docstring.split()[0]
            docstring = ' '.join(docstring.split()[1:])
    #            print "sho:",shortform
    #            print "doc: \"",docstring,"\""
            if docstring != "":
                if docstring[0] == "#":
                    docstring = ""
        ############### END DOCSTRING SHORTFORM ############
        ############### END DOCSTRING SHORTFORM ############

        ############### BEGIN  arguments ###########
        ############### BEGIN  arguments ############
        one = eval(callingscript + "." + deff)
    #        print "one:",one
        two = inspect.getargspec(one)
    #        print 'two:',two
        argumentsin, defaultsin = two.args, two.defaults
    #        print "||| ",argumentsin,defaultsin

        argumentsout = []
        if argumentsin == {}:
            argumentsout = ""
        if argumentsin == []:
            argumentsout = ""
        if argumentsout == []:
            #print ">>||| ",argumentsin," > ",defaultsin
            arguments = zip(argumentsin, defaultsin)
            #print "zip| ",arguments
            arguments = map(list, arguments)
            #print "|||>",arguments

            for line in arguments:
                lineout = str(line[0]) + "=" + str(line[1])
    #                print "|| lineout:",lineout
                argumentsout.append(lineout)

    #        print "$$arguments:",argumentsout[0:],type(argumentsout)
        argumentsout = str(argumentsout)
    #        print "$$ARGUMENTS:",argumentsout,type(argumentsout)
        ############### END  arguments ############
        ############### END  arguments ############

        out.append([deff, shortform, docstring, argumentsout])
    #        print "out:",out
    #    print "OUT:",out
    #    print "  "
    #    print " "

    defs = [i[0] for i in out]
    shortforms = [i[1] for i in out]
    #docstrings = [i[2] for i in out]
    #    print "OUT2:",[defs,shortforms,docstrings]
    if get_duplicate_items(defs) != []:
        print "definition", get_duplicate_items(defs), "is defined more than once"

    #    print "GDI:",get_duplicate_items(shortforms)
    #    print "all shoerforms:", shortforms
    shortforms.append('')
    shortforms.append('')

    if get_duplicate_items(shortforms) != ['']:
        print "shortform for parser:", get_duplicate_items(shortforms), ", is defined more than once"
        print "all shoerforms:", shortforms
        exit()
    return out


def get_duplicate_items(L):
    new = set()
    for item in L:
        if L.count(item) > 1:
            new.add(item)
    return list(new)


def pwd():
    import os
    #return os.getcwd()
    return os.getenv('PWD')


def pwdsplit():
    return pwd().split("/")   # ['', 'nas', 'glensk', 'v', 'PAW_PBE', 'Ag', 'ti_fcc4_bulk', 'low_2x2x2sc_250eV_03x03x03kp_EDIFF1E-1__high_550eV_08x08x08kp']


def pathsplit(path=None):
    if path is None:
        exit("pathsplit needs path to split")
    return str(path).split("/")


def iroundup(x):
    x = x + 0.5
    """iround(number) -> integer
    Round a number to the next integer. 190.0 -> 190; 191.6 -> 192"""
    return int(round(x))


def isdir(path):
    import os
    if os.path.isdir(path) is True:
        return True
    else:
        return False


def isfile(path):
    import os
    if os.path.isfile(path) is True:
        return True
    else:
        return False


def checkdir(path, create=None):
    import os
    import sys
    if os.path.isdir(path) is True:
        return path
    elif create is True:
        os.makedirs(path)
        return path
    else:                         # the case if path does not exist
        import inspect
        if inspect.stack()[1][3] == '<module>':
            print "Could not find path:", path
            print "Check file ", inspect.stack()[1][1], "for checkpath in the main code, not in some definition"
        else:
            print "Could not find path:", path
            print "Check def ", inspect.stack()[1][3], " in file ", inspect.stack()[1][1]
        sys.exit()


def checkfile(path):
    import os
    if os.path.isfile(path) is True:
        return path
    else:
        import inspect
        if inspect.stack()[1][3] == '<module>':
            print "Check file ", inspect.stack()[1][1], "for checkpath in the main code, not in some definition"
        else:
    #            path = nspect.stack()[1][1]
            exit("File \"" + path + "\" does not exist; Check def " +
                 inspect.stack()[1][3] + " in file " + inspect.stack()[1][1])


def exit(error=None):
    import sys
    import inspect
   # from termcolor import colored
    if inspect.stack()[1][3] == '<module>':
        print "Check file ", inspect.stack()[1][1], "for checkpath in the main code, not in some definition"
    else:
        text = ""
        path_full = inspect.stack()[1][1]
        path_name = inspect.getmoduleinfo(path_full).name
        ERROR_GENERAL = "ERROR in Module: " + inspect.stack()[1][3] + "  Skript: " + path_name
        if error is None:
            ERROR = ERROR_GENERAL
        else:
            ERROR = error + " (" + ERROR_GENERAL + ")"
    #            print "KA:",ERROR.split()[0][0:5]
            if ERROR.split()[0][0:5] == "ERROR":
                #text = colored(ERROR, 'red', attrs=['bold'])
                text = printred("ERRORB")
            else:
                #text = colored(ERROR, 'blue', attrs=['bold'])
                #text = printred("ERRORC")
                text = printred(str(error))
        sys.exit(text)
    sys.exit()


def error_write(errortext=None):
    #    print "EE:",errortext
    definition = which_function_calledme()
    file = which_scriptname_isit_whew_i_was_called_from()
    if errortext is None:
    #        return "ERROR def: "+definition+" File: "+file
        return "ERROR> " + file + ".py : " + definition + " : \"" + errortext + "\"" + " <ERROR"

    if errortext is not None:  # es gibt einen errortext
        ## does it contain ERROR> ...<ERROR?
        check = errortext.split("ERROR> ")[-1].split(" <ERROR")[0]
    #        print "check:",check
        if check == errortext:
            return "ERROR> " + file + ".py : " + definition + " : \"" + errortext + "\"" + " <ERROR"
        else:
            return "ERROR> " + check + " <ERROR"
    #        check = errortext.split("ERROR> ")[-1].split(" <ERROR")[0]
    #        print 'check:',check
    #        if check == "ERROR":
    #            out = errortext.split("ERROR")
    #            print "out:",out
    #
    #            return errortext #"LANGERROR: "+file+".py : "+definition+" : \""+errortext+"\""
    #        else:
    #            return "ERROR> "+file+".py : "+definition+" : \""+errortext+"\""+" <ERROR"


def readfile(path=None, line=None):
    if path is None:
        print "i was called from: " + str(which_function_calledme())
        exit(error="Please specify path")

    checkfile(path)
    f = open(path, "r")
    try:
        out = f.read()
    finally:
        f.close()
        if line is None:
            return out.splitlines()
        else:
            return out.splitlines()[line - 1]


def pwdlist():
    return os.getcwd().split("/")   # ['', 'nas', 'glensk', 'v', 'PAW_PBE', 'Ag', 'ti_fcc4_bulk', 'low_2x2x2sc_250eV_03x03x03kp_EDIFF1E-1__high_550eV_08x08x08kp']


    #import locale
    #locale.setlocale(locale.LC_NUMERIC, "")

def format_num(num):
    """Format a number according to given places.
    Adds commas, etc. Will truncate floats into ints!"""
    import locale
    locale.setlocale(locale.LC_NUMERIC, "")
    try:
        inum = int(num)
        return locale.format("%.*f", (0, inum), True)

    except (ValueError, TypeError):
        return str(num)


def get_max_width(table, index):
    """Get the maximum width of the given column index"""
    return max([len(format_num(row[index])) for row in table])


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def is_true_int(s):
    return isinstance(s, int)


def is_string(s):
    return isinstance(s, str)


def is_true_float(s):
    return isinstance(s, float)


def is_list(s):
    return isinstance(s, list)


def is_number(s):
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


def convert_to_number(string):
    if is_int(string) == True:
        value = int(string)
    elif is_float(string) == True:
        value = float(string)
    else:
        print "seems value", string, "is not a number"
        my.exit()
    return value


def parseout(dict):
    def testall():
        pass  # needs to be defined for the loop
    for key in dict:
    #print repr(dict[key]).ljust(10),key.ljust(7)down vote accepted
        if dict[key] is not None and dict[key] != False:
        #            import inspect
        #            print inspect.getargspec(eval(key))
        #            print "key      ---> ",key
        #            print "dict[key] --> ", dict[key]##, len(dict[key])
            if dict[key] == True:
                out = eval(key)()
            elif len(dict[key]) == 0:
                out = eval(key)()
            else:
                out = eval(key)(*dict[key])

            if out is None:
                pass
            else:
                print out

                import inspect
    # functions


def which_function_ami():
    import inspect
    return inspect.stack()[1][3]


def which_function_calledme():
    import inspect
    return inspect.stack()[2][3]


    ############################################################################
    #def which_script_calledme():
    #    import inspect
    #    path = inspect.stack() #[2][1]
    #    print "path:",path
    #    print "path[0]",path[0]
    #    print "path[1]",path[1]
    #    return inspect.getmoduleinfo(path).name
def which_script_show_allinspect():
    import inspect
    return inspect.stack()


def which_scriptpath_isit_whew_i_was_called_from():
    import inspect
    return inspect.stack()[-1][1]


def which_scriptname_isit_whew_i_was_called_from():
    import inspect
    #    return inspect.stack()
    path = inspect.stack()[-1][1]
    #    print ""
    #    print "len:",len(path)
    #    print "-->>>",path
    #    print ""
    return inspect.getmoduleinfo(path).name


def which_scriptpath_isit_whew_i_am():
    import inspect
    #    print "inspect.stack():",inspect.stack()
    path = inspect.stack()[1][1]
    #    print "path",path
    return path


def which_scriptname_isit_whew_i_am():
    import inspect
    #    print "inspect.stack():",inspect.stack()
    path = inspect.stack()[1][1]
    #    print "path",path
    return inspect.getmoduleinfo(path).name
    ############################################################################


def ListFolders(YOUR_PATH):
    import os
    for path in os.listdir(YOUR_PATH):
        if not os.path.isfile(os.path.join(YOUR_PATH, path)):
            yield YOUR_PATH + "/" + path


def numpyarray_to_screenoutput(numpyarray=None):
    if numpyarray is None:
        exit(error="please provide a numpyarray")
    #print numpyarray
    out = numpyarray
    allout = ""
    for line in out:
        #print line,len(line)
        add = "%10.10f  %10.10f  %10.10f" % (line[0], line[1], line[2]) + "\n"
        allout = allout + add
    return allout


def schnittmenge(list1, list2):
    """
    nur elemente die sowolh in list1 als auch in list2 sine
    """
    import myclasses
    sl1 = myclasses.OrderedSetn(list1)
    sl2 = myclasses.OrderedSetn(list2)
    #    print sl1
    #    print sl2
    return list(sl1 & sl2)


def symmetrischdifferenz(list1, list2):
    """
    alle Objekte, die entweder in list1 oder in list2 vorkommen, nicht aber in beiden.
    Anmerkung: set changes order of the lists, OrderedSet nicht!
    """
    import myclasses
    sl1 = myclasses.OrderedSetn(list1)
    sl2 = myclasses.OrderedSetn(list2)
    return list(sl1 ^ sl2)


def differenzmenge(list1, list2):
    """
    Alle Elemente von list1, die nicht in list2 sind
    """
    import myclasses
    sl1 = myclasses.OrderedSetn(list1)
    sl2 = myclasses.OrderedSetn(list2)
    return list(sl1 - sl2)


def vereinigungsmenge(list1, list2):
    """
    vereinige alle Elemente aus list1 und list2
    """
    import myclasses
    sl1 = myclasses.OrderedSetn(list1)
    sl2 = myclasses.OrderedSetn(list2)
    return list(sl1 | sl2)

