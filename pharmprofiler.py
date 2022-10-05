#!/usr/bin/python
# -*- coding: utf-8 -*-

from optparse import OptionParser
from lib.data import *
from lib.patient import *
from lib.reportmaker import *

# --------------------------------------------------------------------------


def main():
    '''
    This main function takes care of the command line input. It is your main hub for generating the pharmacogenetics databases,
    importing patient datasets and report generation.
    Please see the options below and readme for your input options.
    '''
    options = [

        {"short": "n", "long": "name", "action": "store", "dest": "dbname", "default": "database",
         "help": "Specify database name. This will be used to create/access a database. Default is 'database'"},

        {"short": "b", "long": "build", "action": "append", "dest": "dbtasks", "default": [],
         "help":"Download database (specify tables, default is do whole database.\
            Options are pairs (gene-drug pairs), genes, vars (allvars/rsvars/hapvars) \
            and drugs(chemicals"},

        {"short": "a", "long": "action", "action": "append", "dest": "pattasks", "default": [],
         "help":"Which actions to perform on patient file. Current options: 'import' \
        (import important snp sites from vcf), 'haplotype' (phylogenetically haplotype patient),\
        'annotate' (annotate found changes and identified patient haplotypes). Use 'guidelines' to find patient matching guidelines and add to DB.\
        Use 'annotate' to do the same for annotations. Use 'report' after this to generate a patient report using this info (use -o to specify a location for the resulting .tex file.)"},

        {"short": "p", "long": "patient", "action": "store", "dest": "gvcf", "default": None,
         "help": "Patient compressed vcf [g.vcf.gz] file to parse"},

        {"short": "o", "long": "outfile", "action": "store", "dest": "outfile", "default": "",
         "help": "Location to send report TEX file to."},

        {"short": "r", "long": "reset", "action": "append", "dest": "reset", "default": [],
         "help":"Database table to reset (mostly for testing purposes)"}]

    # Create object that parses command line options
    parser = OptionParser(usage='usage: %prog [options] filename',
                          version='%prog 1.0')
    for o in options:
        parser.add_option(
            '-{}'.format(o['short']),
            '--{}'.format(o['long']),
            action=o['action'],
            dest=o['dest'],
            default=o['default'],
            help=o['help']
        )

    # Get command line arguments
    (options, args) = parser.parse_args()

    if len(options.dbtasks) > 0:
        # Check if any database actions have been selected - if so, execute these actions.
        create_db(options.dbname, options.dbtasks, options.reset)
    if options.gvcf:
        # Check if a GVCF has been supplied - if so, create patient with given information
        create_patient(options.dbname, options.gvcf, options.pattasks, options.outfile)


def create_db(dbname, tables, reset):
    '''
    This function takes care of database creation. It executes various functions necessary
    for table creation, PharmGKB API communication and data insertion.
    '''
# --------------------------------------------------------------------------
    # Create database creating object, DataCollector
    d = DataCollector(dbname)
    # Check if any tables need to be resetted (deleted and re-created), mostly for testing purposes
    if len(reset) > 0:
        for table in reset:
            d.remake_table(reset)
# --------------------------------------------------------------------------
    # This dictionary stores the functions coupled to the command line options

    options = OrderedDict([
        ("pairs", d.get_pairs),
        ("drugs", d.get_drugs),
        ("drugmatches", d.get_drug_matches),
        ("drugvars", d.get_drug_vars),
        ("genes", d.get_gene_data),
        ("haplotypes", d.get_haplotypes),
        ("hapvars", d.get_hap_vars),
        ("etcvars", d.get_non_rs),
        ("indels", d.convert_indels),
        ("annotations", d.get_annotations),
        ("guidelines", d.get_guidelines),
        ("guideoptions", d.get_guide_options),
        ("bed", d.bed_file)
    ])

    # Create options that do multiple things in a row (for example creating all variant-related tables)
    options['vars'] = [options['drugvars'], options["hapvars"], options["etcvars"]]
    options['all'] = options.values()
# --------------------------------------------------------------------------
    # Go through tables supplied in command line
    for table in tables:
        # Get selected tables
        o = options[table]
        '''
        Check if an API authentication token is present
        (needed for getting variant annotations from PharmGKB API)
        otherwise create one and add to the DataCollector object.
        '''
        if not hasattr(d, "authobj"):
            d.authenticate()
        # Execute functions matching your selected options
        if type(o) is not list:
            options[table]()
        elif type(o) is list:
            for item in o:
                item()
    # Commit changes to database
    d.conn.commit()


def create_patient(dbname, gvcf, tables, outfile):
    # Create patient object
    p = Patient(dbname, gvcf)
    # Create interpreter object
    i = Interpreter(p)
    # Create report making object
    r = ReportMaker(p, outfile)

    # Dictionary with functions related to patient data parsing
    options = OrderedDict([
        ("import", p.get_positions),
        ("haplotype", p.get_haplotypes),
        ("score", p.hap_scorer),
        ("genotype", i.get_genotype),
        ("guidelines", i.find_guidelines),
        ("annotate", i.make_annotation),
        ("report", [r.make_json, r.make_report])
    ])

    # Create options that do multiple things in a row
    options['all'] = options.values()
# ----------------------------------------------
    # Go through patient tables that need to be created
    # and execute matching function from above dictionary
    for table in tables:
        try:
            o = options[table]
            if type(o) is not list:
                options[table]()
            elif type(o) is list:
                for item in o:
                    item()
        except:
            raise
            print "Invalid option entered. \n Valid options: {}".format(", ".join(options.keys()))

    # Commit changes to database
    p.conn.commit()

# --------------------------------------------------------------------------------


# Make sure this only runs when you actively execute it
# not when sourcing this file
if __name__ == '__main__':
    main()
