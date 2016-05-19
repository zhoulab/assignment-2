import re
import os
import logging
import time
import csv

from bs4 import BeautifulSoup
from goatools import obo_parser

# Regular expresions for terms relating to cell death
RE_CELL_DEATH = [re.compile('death', re.I), re.compile('apopto', re.I),
                 re.compile('caspase', re.I)]

# We are filtering out GO "levels" higher than this value (inclusive)
LEVEL = 3

ALPHA = 1e-3

OBO_FILE = 'go-basic.obo'
OUTPUT_FILE = "go_results.txt"


def get_rows(soup):
    """
        Parameters
        ----------
        soup: BeautifulSoup object
            Contains a parsed 'geneOntology.html' file

        Returns
        -------
        ({head1: val1-1, head2: val1-2, ...},
         {head1: val2-1, head2: val2-2, ...}, ...)
                head(j) is the jth header and val(i)-(j) is
                the value under head(j) for the i-th row
    """
    head = [col.text for col in soup.find('tr').findChildren()]
    row_values = [[col.text for col in row.findChildren()]
                  for row in (soup.findAll('tr')[1:])]
    for row in row_values:
        yield dict(zip(head, row))


def print_header():
    """Prints a formatted header to match row output"""
    print '%-21s%-43s%-11s%-20s%-12s%s' % ('Filename', 'Term',
                                           'P-value', 'GO Tree',
                                           '# of Genes',
                                           '# of Target Genes')


def print_row(row):
    """Prints only select columns from the row in a formatted output"""
    print '%-21s%-43s%-11s%-20s%-12s%s' % (row['Filename'], row['Term'],
                                           row['P-value'], row['GO Tree'],
                                           row['# of Genes in Term'],
                                           row['# of Target Genes in Term'])


def get_ancestors(go_obj):
    """
        Returns
        -------
        {ancestor1, ancestor2, ancestor3, ...} : set
    """

    if not go_obj.parents:
        return set()

    return set(go_obj.parents) | set.union(*[get_ancestors(parent)
                                             for parent in go_obj.parents])


if __name__ == "__main__":
    """
        Goes through each folder in the 'GO_folders' directory
        and analyzes the 'geneOntology.html' file.
        'GO_folders' directory must be a sibling of this file's parent directory
    """
    log = logging.getLogger()
    handler = logging.StreamHandler()
    log.addHandler(handler)
    log.setLevel(logging.DEBUG)

    log.info('Parsing %s for level checking', OBO_FILE)
    p = obo_parser.GODag(OBO_FILE)

    os.chdir('..')
    os.chdir('DmGOs')
    folders = next(os.walk('.'))[1][1:3]

    start = time.clock()

    rows = []
    for i, folder in enumerate(folders):
        log.info('Searching in folder %i of %i: %s',
                 i + 1, len(folders), folder)
        os.chdir(folder)
        with open('geneOntology.html') as go_file:
            soup = BeautifulSoup(go_file, "html.parser")
            for row in get_rows(soup):
                # filter using level, pvalue. some GO IDs are not valid
                # (ex. 'fly-chr-') so we will ignore those
                if ('GO:' in row['GO ID'] and p[row['GO ID']].level > LEVEL and
                        row['P-value'] > ALPHA):
                    GOTerm_obj = p[row['GO ID']]
                    row['Filename'] = folder
                    row['Level'] = GOTerm_obj.level
                    row['Depth'] = GOTerm_obj.depth
                    row['Ancestors'] = [a.name for a in get_ancestors(GOTerm_obj) if a.level > LEVEL]
                    rows.append(row)
            log.info('Found %i significant terms in %s', len(rows), folder)
        os.chdir('..')

    rows = sorted(rows, key=lambda k: float(k['P-value']))
    print 'Found %i terms with p-value > %f:' % (len(rows), ALPHA)

    os.chdir('..')
    with open(OUTPUT_FILE, "w") as out_file:
        print 'Writing file...'
        filewriter = csv.writer(out_file, delimiter='\t')
        header = ['Filename', 'Term', 'Level', 'Depth', 'P-value',
                  '# of Genes in Term', '# of Target Genes in Term',
                  '# of Total Genes', '# of Target Genes', 'Ancestors']
        filewriter.writerow(header)
        for row in rows:
            filewriter.writerow([row[head] for head in header])

    end = time.clock()
    print 'Run time was %4.2fs' % (end - start)
