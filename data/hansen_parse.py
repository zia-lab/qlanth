#/usr/bin/env python3
import pandas as pd
from collections import Counter
from pandas import ExcelWriter
import os
import numpy as np

info = '''┌────────────────────────────────────────────────────────────────────────┐
│~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~│
│~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~│
│~~~~┌──────────────────────────────────────────────────────────────┐~~~~│
│~~~~│  This Python script parses and cleans the raw data from the  │~~~~│
│~~~~│                          tables of                           │~~~~│
│~~~~│ Hansen, JE, BR Judd, and Hannah Crosswhite. “Matrix Elements │~~~~│
│~~~~│ of Scalar Three-Electron Operators for the Atomic f-Shell.”  │~~~~│
│~~~~│ Atomic Data and Nuclear Data Tables 62, no. 1 (1996): 1–49.  │~~~~│
│~~~~│ It creates an excel spreadsheet where each sheet corresponds │~~~~│
│~~~~│  to an f^N configuration, and it also produces a single csv  │~~~~│
│~~~~│    file with a column that allows discriminating for the     │~~~~│
│~~~~│                  different configurations.                   │~~~~│
│~~~~│        The abstract from that paper is the following:        │~~~~│
│~~~~│                                                              │~~~~│
│~~~~└──────────────────────────────────────────────────────────────┘~~~~│
│~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~│
│~~~~~~~~                                                        ~~~~~~~~│
│~~~~~~~~   "Tables are provided for the matrix elements of an   ~~~~~~~~│
│~~~~~~~~  orthogonal set of Hermitian three-electron operators  ~~~~~~~~│
│~~~~~~~~  ti for the states of the f shell. The ti are scalar   ~~~~~~~~│
│~~~~~~~~   with respect to the total spin S and total orbital   ~~~~~~~~│
│~~~~~~~~  angular momentum L, and they are among the effective  ~~~~~~~~│
│~~~~~~~~    operators needed to be included in an f-electron    ~~~~~~~~│
│~~~~~~~~ Hamiltonian in order to represent the coupling of the  ~~~~~~~~│
│~~~~~~~~ ground configuration f N to excited configurations via ~~~~~~~~│
│~~~~~~~~       the inter-electronic Coulomb interaction."       ~~~~~~~~│
│~~~~~~~~                                                        ~~~~~~~~│
│~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~│
│~~~~┌──────────────────────────────────────────────────────────────┐~~~~│
│~~~~│                          IMPORTANT:                          │~~~~│
│~~~~│                                                              │~~~~│
│~~~~│  + The values used here for t_2 are instead those of t'_2.   │~~~~│
│~~~~│   + The .xls file was produced with help of OCR, but their   │~~~~│
│~~~~│  values were later tested against another version of these   │~~~~│
│~~~~│           tables and no discrepancies were found.            │~~~~│
│~~~~└──────────────────────────────────────────────────────────────┘~~~~│
│~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~│
│~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~│
└────────────────────────────────────────────────────────────────────────┘'''

tindices = [2,3,4,6,7,8,11,12,14,15,16,17,18,19]
tsymbols = ['t%d' % s for s in tindices]
workdir = os.getcwd()
export_fname = os.path.join(workdir, 'ti-hansen1996.xls')
raw_excel = os.path.join(workdir, 'hansen.xlsx')
make_plots = True


def num_cleaner(thing):
    thing = str(thing)
    thing = thing.replace(' ','')
    thing = thing.replace('이','01')
    return float(thing)

def term_cleaner(numE, thing):
    ambiguities = '51 5I-61 6I-20 2O-30 3O-IP 1P-IF 1F-IN 1N-IQ 1Q'
    ambiguities = [tuple(s.split(' ')) for s in ambiguities.split('-')]
    thing = str(thing)
    thing = thing.replace('工','I')
    if numE == 7:
        thing = thing.replace('2F0','2F10')
        thing = thing.replace('2G0','2G10')
    for k in (1,2,3,4):
        for i in range(k*100+10,k*100+20):
            thing = thing.replace(str(i),'%dI%s' % (k,str(i)[2]))
    for (l, r) in ambiguities:
        thing = thing.replace(l, r)
    return thing

def indices(lst, item):
    return [i for i, x in enumerate(lst) if x == item]

def parse_hansen(saving=True, verbose=False, make_plots=False):
    hansen = {}
    for numElectrons in range(3,8):
        # read the sheet in the excel spreadsheet
        hanFrame = pd.read_excel(raw_excel, sheet_name=numElectrons-3)
        # tidy up the term symbols
        hanFrame["bterm"] = hanFrame["bterm"].apply(lambda x: term_cleaner(7,x))
        hanFrame["kterm"] = hanFrame["kterm"].apply(lambda x: term_cleaner(7,x))
        # tidy up the numbers
        for tindex in tindices:
            hanFrame['t%d' % tindex] = hanFrame['t%d' % tindex].apply(num_cleaner)
        hansen[numElectrons] = hanFrame

    # There are some repeated rows in the tables by Hansen 
    # this loop is related to that
    for numElectrons in range(3,8):
        hanFrame = hansen[numElectrons]
        termpairs = list(zip(hanFrame['bterm'], hanFrame['kterm']))
        if verbose:
            print(len(termpairs), len(set(termpairs)))
        counts = Counter(termpairs)
        for k,v in counts.items():
            if v == 1:
                continue
            reps = indices(termpairs, k)
            reprows = []
            for rep in reps:
                reprows.append(tuple(hanFrame.iloc[rep]))
            if verbose:
                print(len(reprows),len(set(reprows)))

    # Remove those duplicate rows
    for numElectrons in range(3,8):
        hanFrame = hansen[numElectrons]
        newFrame = hanFrame.drop_duplicates()
        hansen[numElectrons] = newFrame

    # add the symmetric elements
    for numElectrons in range(3,8):
        hanFrame = hansen[numElectrons]
        termpairs = list(zip(hanFrame['bterm'], hanFrame['kterm']))
        nondiag = [tp[0] != tp[1] for tp in termpairs]
        diag = [tp[0] == tp[1] for tp in termpairs]
        nonDiagFrame = hanFrame[nondiag]
        diagFrame = hanFrame[diag]
        extraFrame = nonDiagFrame.rename(columns={'bterm':'kterm','kterm':'bterm'})
        hansen[numElectrons] = pd.concat([nonDiagFrame, diagFrame, extraFrame])
        hansen[numElectrons].sort_values('bterm', inplace=True)

    def save_xls(list_dfs, xls_path):
        with ExcelWriter(xls_path) as writer:
            for n, df in enumerate(list_dfs):
                if n > 0:
                    df.to_excel(writer,'N = %s' % (n+2), index=False)
                else:
                    df.to_excel(writer,'comments',  header=False, index=False)

    comment = 'This is a digital version of the tables from Hansen, JE, BR Judd, and Hannah Crosswhite. “Matrix Elements of Scalar Three-Electron Operators for the Atomic f-Shell.” Atomic Data and Nuclear Data Tables 62, no. 1 (1996): 1–49. The were digitized via OCR, take them with a grain of salt.'
    export = pd.DataFrame([[comment]])
    export = list(hansen.values())
    export = [pd.DataFrame([[comment]])] + export
    if saving:
        print("Saving to %s ..." % export_fname)
        save_xls(export, export_fname)

    for numElectrons in range(3,8):
        hanFrame = hansen[numElectrons]
        hanFrame['f^N'] = [numElectrons]*len(hanFrame)
        hansen[numElectrons] = hanFrame

    bigHansenFrame = pd.concat(hansen)
    bigHansenFrame = bigHansenFrame[['f^N','bterm', 'kterm', 't2', 't3', 't4', 't6', 't7', 't8', 't11', 't12',
        't14', 't15', 't16', 't17', 't18', 't19']]
    csv_fname = os.path.join(workdir, 'ti-hansen1996.csv')
    if saving:
        print("Saving to %s ..." % csv_fname)
        bigHansenFrame.to_csv(csv_fname, index=False)
    if make_plots:
        import cmasher as cmr
        from matplotlib import pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        tstrings = 't₂,t₃,t₄,t₆,t₇,t₈,t₁₁,t₁₂,t₁₄,t₁₅,t₁₆,t₁₇,t₁₈,t₁₉'.split(',')
        with PdfPages('ti-figs.pdf') as pdf:
            for numElectrons in [3,4,5,6,7]:
                print("Writing plots for f^%d:" % numElectrons)
                for tsindex, ts in enumerate(tsymbols):
                    hanFrame = bigHansenFrame[bigHansenFrame['f^N']==numElectrons]
                    terms = list(hanFrame['bterm'].unique())
                    hanDict = hanFrame.set_index(['bterm', 'kterm']).to_dict()
                    figwidth = (25/4*(numElectrons-3) + 5)
                    x_labels = terms
                    y_labels = terms
                    fonsize = (60/4*(numElectrons-3) + 20)
                    tiMatrix = [[hanDict[ts].get((bterm, kterm),0.) for bterm in terms] for kterm in terms]
                    tiMatrix = np.array(tiMatrix)
                    plt.figure()
                    minv, maxv = np.max(tiMatrix), np.min(tiMatrix)
                    valuerange = max(abs(minv), abs(maxv))
                    fig, ax = plt.subplots(figsize=(figwidth, figwidth))
                    ax.imshow(tiMatrix, 
                            vmin = -valuerange,
                            vmax = valuerange, 
                            cmap = cmr.watermelon)
                    # Set the ticks - the positions where the labels should appear
                    ax.set_xticks(np.arange(len(x_labels)))
                    ax.set_yticks(np.arange(len(y_labels)))
                    
                    # Set the labels and rotate x labels
                    ax.set_xticklabels(x_labels, rotation=90)
                    ax.set_yticklabels(y_labels)
                    
                    # Add ticks to the bottom and left Axis as well
                    ax.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True, labeltop=True)
                    ax.tick_params(axis='y', which='both', right=True, left=True, labelright=True, labelleft=True)
                    ax.set_title('~%s in f%s~' % (tstrings[tsindex], '³⁴⁵⁶⁷'[numElectrons-3]),fontsize=fonsize)
                    figname = 'f%d-t%d.pdf' % (numElectrons, tindices[tsindex])
                    plt.tight_layout()
                    pdf.savefig()
                    plt.close()
    return bigHansenFrame

if __name__ == '__main__':
    print(info)
    parse_hansen(saving=True, make_plots=make_plots)