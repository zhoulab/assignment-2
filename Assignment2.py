import re
import os
import logging
import time

from bs4 import BeautifulSoup
from goatools import obo_parser

# Regular expresions for terms relating to cell death
RE_CELL_DEATH = [re.compile('death', re.I), re.compile('apopto', re.I),
                 re.compile('caspase', re.I)]

# We are filtering out GO "levels" higher than this value (inclusive)
LEVEL = 3

OBO_FILE = 'go-basic.obo'


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
    folders = next(os.walk('.'))[1]

    cell_death_rows = []
    significant_rows = []

    start = time.clock()

    for i, folder in enumerate(folders):
        log.info('Searching in folder %i of %i: %s',
                 i + 1, len(folders), folder)
        os.chdir(folder)
        with open('geneOntology.html') as f:
            soup = BeautifulSoup(f, "html.parser")
            # count cell death rows in each file for log info
            count = 0
            for row in get_rows(soup):
                # filter based on level. some GO IDs are not actually GO IDs
                # (ex. 'fly-chr-') so we will ignore those
                if 'GO:' in row['GO ID'] and p[row['GO ID']].level > LEVEL:
                    row['Filename'] = folder
                    if any(rgx.match(row['Term']) for rgx in RE_CELL_DEATH):
                        cell_death_rows.append(row)
                        count += 1
                    if row['P-value'] < 1e-15:
                        significant_rows.append(row)
            log.info('Found %i cell death related terms in %s', count, folder)
        os.chdir('..')

    cell_death_rows = sorted(cell_death_rows,
                             key=lambda k: float(k['P-value']))
    print 'Found %i terms relating to cell death:' % len(cell_death_rows)
    print_header()
    for row in cell_death_rows:
        print_row(row)
    print ''
    print 'Found %i terms with p-value < 1e-15:' % len(significant_rows)
    for row in significant_rows:
        print_row(row)

    end = time.clock()
    print 'Run time was %4.2fs' % (end - start)
