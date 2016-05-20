import re
import os
import logging
import time
import csv
import sys

from bs4 import BeautifulSoup
from goatools import obo_parser

# Regular expresions for terms relating to cell death
RE_CELL_DEATH = [re.compile('death', re.I), re.compile('apopto', re.I),
                 re.compile('caspase', re.I)]

# We are filtering out GO "levels" higher than this value (inclusive)
LEVEL = 3

ALPHA = 1e-3
TOP_RESULTS_CAP = 10

OBO_FILE = 'go-basic.obo'

HEADER = ['GO ID', 'Term', 'Level', 'Depth', 'P-value',
          '# of Genes in Term', '# of Target Genes in Term',
          '# of Total Genes', '# of Target Genes', 'Common Genes']

log = logging.getLogger()
handler = logging.StreamHandler()
log.addHandler(handler)
log.setLevel(logging.DEBUG)


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


def get_ancestors(go_obj, i=1):
    """
        Start recursion with i=1 to get desired diffs

        Returns
        -------
        {(ancestor1, diff), (ancestor2, diff), ...} : set of 2-tuples
    """

    if not go_obj.parents:
        return set()

    return set([(parent, i)
                for parent in go_obj.parents]) | set.union(*[get_ancestors(parent, i + 1) for parent in go_obj.parents])


def generate_parent_levels(rows):
    """
        Parameters
        ----------
        rows: list of dicts
            GO rows containing important info (see HEADER)

        Returns
        -------
        int
            Highest level difference between an ancestor and its term
    """
    max_diff = -1
    for row in rows:
        max_diff_in_row = max([diff for (__, diff) in row['Ancestors']])
        if max_diff_in_row > max_diff:
            max_diff = max_diff_in_row
        for i in range(1, (max_diff_in_row + 1)):
            log.debug('creating column for "Parent %i" in %s', i, row['Term'])
            row['Parents (diff: ' + str(i) + ")"] = []
        for ancestor, diff in row['Ancestors']:
            row['Parents (diff: ' + str(diff) + ")"].append(ancestor.name)
    return max_diff


def get_updated_header(header, max_diff):
    """Takes in a header and returns an updated list based on max difference"""
    new_header = header
    for i in range(1, max_diff + 1):
        new_header.append('Parents (diff: ' + str(i) + ")")
    return new_header


def memory_usage_resource():
    import resource
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem


def assignment2():
    """
        Goes through each folder in the 'GO_folders' directory
        and analyzes the 'geneOntology.html' file.
        'GO_folders' directory must be a sibling of this file's parent directory
    """
    log.debug('Memory usage: %s', memory_usage_resource())
    log.info('Parsing %s for level checking', OBO_FILE)
    p = obo_parser.GODag(OBO_FILE)
    log.debug('Memory usage: %s', memory_usage_resource())

    os.chdir('..')
    os.chdir('DmGOs')
    folders = next(os.walk('.'))[1]

    start = time.clock()

    for i, folder in enumerate(folders):
        log.info('Searching in folder %i of %i: %s',
                 i + 1, len(folders), folder)
        log.debug('Memory usage: %s', memory_usage_resource())
        os.chdir(folder)

        rows = []
        log.debug('Opening geneOntology.html')
        with open('geneOntology.html') as go_file:
            log.debug('Memory usage: %s', memory_usage_resource())
            log.debug('Parsing geneOntology.html with BeautifulSoup')
            soup = BeautifulSoup(go_file, "html.parser")
            time.sleep(1)
            log.debug('Memory usage: %s', memory_usage_resource())
            for row in get_rows(soup):
                # filter using level, pvalue. some GO IDs are not valid
                # (ex. 'fly-chr-') so we will ignore those
                if ('GO:' in row['GO ID'] and p[row['GO ID']].level > LEVEL and
                        row['P-value'] > ALPHA):
                    GOTerm_obj = p[row['GO ID']]
                    row['Filename'] = folder
                    row['Level'] = GOTerm_obj.level
                    row['Depth'] = GOTerm_obj.depth
                    row['Ancestors'] = [(a, d) for (a, d) in get_ancestors(GOTerm_obj) if a.level > 2]
                    rows.append(row)
            soup = None
            log.info('Found %i significant terms (p<%f) in %s', len(rows), ALPHA, folder)
            rows = sorted(rows, key=lambda k: float(k['P-value']))[:TOP_RESULTS_CAP]
            max_diff = generate_parent_levels(rows)
            new_header = get_updated_header(HEADER, max_diff)
            # make sure each row has all keys, regardless if needed.
            for row in rows:
                for key in new_header:
                    if key not in row:
                        row[key] = None

            # change directory to parent of DmGOs
            os.chdir('..')
            os.chdir('..')
            with open(folder + '-results.txt', 'w') as out_file:
                filewriter = csv.writer(out_file, delimiter='\t')
                log.info('Writing output file for %s (top %i significant terms)', folder, TOP_RESULTS_CAP)
                filewriter.writerow(new_header)
                for row in rows:
                    for col in new_header:
                        if 'Parents' in col and row[col]:
                            row[col] = ', '.join(row[col])
                    filewriter.writerow([row[col] for col in new_header])
            log.debug('Memory usage: %s', memory_usage_resource())
            os.chdir('DmGOs')

    end = time.clock()
    print 'Run time was %4.2fs' % (end - start)

if __name__ == "__main__":
    assignment2()
