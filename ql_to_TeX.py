#!/usr/bin/env python3

import os
import textwrap

# parses qlanth.m, and fittings.m into blocks of code

for fname in ['qlanth.m','fittings.m']:
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
        with open(fname,'w') as f:
            f.write(textwrap.dedent(usage_and_def))

