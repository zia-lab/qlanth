#!/usr/bin/env python3

import os
import textwrap
import shutil

# parses qlanth.m, and fittings.m into blocks of code

thesis_fundef_folder = '/Users/juan/david-thesis/fundefs'
thesis_code_folder = '/Users/juan/david-thesis/code'
code_files = ['qlanth.m', 'fittings.m', 'misc.m', 'qalculations.m', 'qonstants.m', 'qplotter.m',
              'fitting-LiYF4-and-LaF3.m']

figures_from_to = [
    ('/Users/juan/ZiaLab/Codebase/qlanth/doc/pr_code.png',
     '/Users/juan/david-thesis/figures/pr_code.png'),
    ('/Users/juan/ZiaLab/Codebase/qlanth/doc/ion_in_lattice.jpg',
    '/Users/juan/david-thesis/figures/ion_in_lattice.jpg'),
    ('/Users/juan/ZiaLab/Codebase/qlanth/doc/nd_code.png',
     '/Users/juan/david-thesis/figures/nd_code.png')
              ]

for (source, dest) in figures_from_to:
    print("Copying %s to %s" % (source, dest))
    shutil.copy(source, dest)

for fname in ['qlanth.m','fittings.m']:
    print(fname)
    data = open(fname,'r').read()
    functions = {}
    lines = data.split('\n')
    idx = 0
    while True:
        line = lines[idx]
        if '::usage'  in line and 'usageTemplate' not in line:
            function = line.split('::usage')[0].strip().split(':')[0].strip()
            usagelines = []
            # parse the usage string
            while True:
                line = lines[idx].strip()
                usagelines.append(lines[idx])
                if lines[idx] == '' or '";' in line:
                    break
                idx += 1
            usage_string = '\n'.join(usagelines)
            deflines = []
            while True:
                idx += 1
                line = lines[idx].strip()
                deflines.append(lines[idx])
                if lines[idx] == '':
                    break
            def_string = '\n'.join(deflines)
            usage_string = usage_string.replace('\\"',"''")
            def_string = def_string.replace('\\"',"''")
            functions[function] = (usage_string, def_string)
        idx += 1
        if idx == len(lines):
            break

    line_template = '''\\begin{lstlisting}
    %s
    %s
    \\end{lstlisting}'''
    for key, (usage, defn) in functions.items():
        usage_and_def = "%s\n%s" % (usage,defn)
        fname = os.path.join('./doc','fundefs',key+'.tex')
        print("Writing to %s" % fname)
        writing_this = textwrap.dedent(usage_and_def)
        with open(fname,'w') as f:
            f.write(writing_this)
        fname2 = os.path.join(thesis_fundef_folder,key+'.tex')
        print("Writing to %s" % fname2)
        with open(fname2,'w') as f:
            f.write(writing_this)

for code_file in code_files:
    # copy to the code folder
    # copy to the thesis code folder
    # copy to the thesis code folder
    source = os.path.join('./',code_file)
    dest = os.path.join(thesis_code_folder,code_file)
    shutil.copy(source, dest)
         

