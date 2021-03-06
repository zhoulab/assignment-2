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
# We are filtering GO terms by level
LEVEL = 4
ALPHA = 1e-3
TOP_RESULTS_CAP = 10
HEADER = ['Term', 'P-value', 'Level_4_Traceback']

BASE_DIR = os.path.dirname(os.getcwd())
OBO_FILE = os.path.join(BASE_DIR, 'data/go-basic.obo')

GO_FOLDERS = [os.path.join(BASE_DIR, 'data/GOs', f)
              for f in ['DmGOs', 'MammalGOs']]
OUTPUT_DIRECTORY = os.path.join(BASE_DIR, 'results/P53-ChIPSeq-GO-results')
if not os.path.exists(OUTPUT_DIRECTORY):
    os.makedirs(OUTPUT_DIRECTORY)


def get_rows(soup):
    """
        Yield rows in a dictionary structure.

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


def get_ancestors(go_obj, i=1):
    """
        Recursively get all ancestors of a GOTerm object.

        Returns
        -------
        {(ancestor1, diff), (ancestor2, diff), ...} : set of 2-tuples

        Notes
        -----
        Starts recursion with i=1 to get desired diffs
    """

    if not go_obj.parents:
        return set()

    return set([(parent, i) for parent in go_obj.parents]) | \
        set.union(*[get_ancestors(parent, i + 1) for parent in go_obj.parents])


def generate_parent_levels(rows):
    """
        Generate parent keys for each row.

        Parameters
        ----------
        rows: list of dicts
            GO rows containing important info (see HEADER)

        Returns
        -------
        int
            Highest level difference between an ancestor and its term.
            Used by get_updated_header()
    """
    for row in rows:
        row['Level_4_Traceback'] = []
        if row['Level'] == 4:
            row['Level_4_Traceback'].append(row['Term'])
        if row['Ancestors']:
            max_diff_in_row = max([diff for (__, diff) in row['Ancestors']])
            for ancestor, diff in row['Ancestors']:
                if diff == max_diff_in_row:
                    row['Level_4_Traceback'].append(ancestor.name)


def get_updated_header(header, max_diff):
    """Takes in a header and returns an updated list based on max difference"""
    return header + ['Parents (diff: ' + str(i) + ")"
                     for i in range(1, max_diff + 1)]


def get_significant_rows(soup, p, alpha, level_cutoff):
    """
        Yield significant rows based on alpha, with additional keys attached.

        Parameters
        ----------
        soup: BeautifulSoup object
            Contains a parsed 'geneOntology.html' file
        p: GODag object
            data from .obo file parsed by obo_parser

        Returns
        -------
        ({row1}, {row2}, ...)
            see get_rows() for row data structure
    """
    for row in get_rows(soup):
        # filter using level, pvalue. some GO IDs are not valid
        # (ex. 'fly-chr-') so we will ignore those
        if ('GO:' in row['GO ID'] and
                p[row['GO ID']].level >= level_cutoff and
                row['P-value'] > alpha):
            term_obj = p[row['GO ID']]
            row['Level'] = term_obj.level
            row['Depth'] = term_obj.depth
            row['Ancestors'] = [(a, d) for (a, d) in get_ancestors(term_obj)
                                if a.level >= level_cutoff]
            yield row


def memory_usage_resource():
    """
        Tool for debugging memory usage.
        Source: http://fa.bianp.net/
    """
    import resource
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem


def P53_ChIPSeq_GO():
    """
        Go through each sample folder in 'DmGOs' and 'MammalGOs'
        and analyze the 'geneOntology.html' file.
    """
    log.info('Parsing %s for level checking', OBO_FILE)
    p = obo_parser.GODag(OBO_FILE, load_obsolete=True)

    start = time.clock()

    for go_folder in GO_FOLDERS:
        sample_folders = next(os.walk(go_folder))[1]

        summary_rows = [[] for i in range(TOP_RESULTS_CAP)]

        for i, sample in enumerate(sample_folders):
            log.info('Searching in sample %i of %i: %s', i + 1, len(sample_folders), sample)
            with open(os.path.join(go_folder, sample, 'geneOntology.html')) as go_file:
                soup = BeautifulSoup(go_file, "html.parser")
                rows = get_significant_rows(soup, p, ALPHA, level_cutoff=LEVEL)
            rows = sorted(list(rows), key=lambda k: float(k['P-value']))[:TOP_RESULTS_CAP]
            generate_parent_levels(rows)

            with open(os.path.join(OUTPUT_DIRECTORY, sample + '-results.txt'), 'w') as out_file:
                filewriter = csv.writer(out_file, delimiter='\t')
                log.info('Writing output file for %s (top %i significant terms)',
                         sample, TOP_RESULTS_CAP)
                filewriter.writerow(HEADER)
                for i, row in enumerate(rows):
                    # make parent lists pretty
                    for col in HEADER:
                        if col == 'Level_4_Traceback':
                            row[col] = ', '.join(row[col])
                    filewriter.writerow([row[col] for col in HEADER])
                    summary_rows[i] += [row[col] for col in HEADER]

        log.info("Writing summary file")
        summary_header = []
        for sample in sample_folders:
            summary_header += [sample + "_" + val for val in HEADER]

        with open(os.path.join(OUTPUT_DIRECTORY, os.path.basename(go_folder) + '-summary.txt'), 'w') as f:
            log.info(f)
            filewriter = csv.writer(f, delimiter='\t')
            filewriter.writerow(summary_header)
            for row in summary_rows:
                filewriter.writerow(row)

    end = time.clock()
    log.info('Run time was %4.2fs', end - start)
    log.info('Memory usage: %s', memory_usage_resource())


if __name__ == "__main__":
    log = logging.getLogger()
    handler = logging.StreamHandler()
    log.addHandler(handler)
    log.setLevel(logging.INFO)

    if not os.path.exists(OUTPUT_DIRECTORY):
        os.makedirs(OUTPUT_DIRECTORY)

    P53_ChIPSeq_GO()
