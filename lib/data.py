#!/usr/bin/python
# -*- coding: utf-8 -*-

from collections import OrderedDict
from modules.pgkb_functions import authenticate, hg19conv, get_json, get_ref, deconstruct_guideline
import urllib2
import json
from tqdm import tqdm
from db import Database
import pprint as pp
import os
import sys
import csv
import sqlite3
import time

# ---------------------------------------------------------------------


class DataCollector(Database):
    '''
    This class collects data from the PharmGKB API
    and uses it to fill the local database that is used later on.
    '''

    def __init__(self, dbname):
        """
        Initiates database object (see 'db.py')
        """
        # Get current location
        self.path = os.path.dirname(__file__)
        if len(self.path) == 0:
            self.path = os.getcwd()

        # Generate folder location and database file name
        dbfolder = os.path.join(self.path, 'db')
        self.dbpath = os.path.join(dbfolder, '%s.db' % dbname)

        # Initiate database object at the above filepath (db.py)
        Database.__init__(self, self.dbpath)

    def authenticate(self):  # This function has the same name as in "pgkb_functions.py"
        # Get token from PharmGKB (REQUIRES ACCOUNT)
        self.authobj = authenticate()

    def get_pairs(self):
        """
        Create table for drug-gene pairs, to be used later to fetch drug data
        :return:
        """
        # Remove and re-create table base (initiate empty table)
        self.remake_table('pairs')
        print 'Getting gene-drug pairs!'
        # Define uri for finding well-annotated pairs (given by pharmgkb)
        uri = 'https://api.pharmgkb.org/v1/report/selectPairs'
        # Make an API call to 'uri' to get list of gene-drug pairs
        data = get_json(uri, self.authobj)
        # Fill Jinja template with recieved data
        # see DB object and template folder for details
        template = self.insert_sql("pairs")

        # For each item in recieved JSON dictionary
        for response in tqdm(data):
            # Fill INSERT SQL template
            sql = template.render(json=response)
            # Execute SQL query
            self.sql.executescript(sql)
            # Commit changes
            self.conn.commit()

    def get_drugs(self):
        '''
        Get list of all available drugs from downloaded TSV file
        and insert into Drugs table.
        '''
        # RECOMMENDED: AUTO-DOWNLOAD NEW DRUGS TSV FILES
        # SEE: https://www.pharmgkb.org/downloads

        # Remove and re-create table base (initiate empty table)
        self.remake_table("drugs")
        # Make sure table doesn't get cut off due to char limit
        csv.field_size_limit(sys.maxint)

        # Open CSV file and go through rows
        with open(self.path+"/drugs/drugs.tsv", "rb") as f:
            # Parse CSV
            tsv = csv.reader(f, delimiter='\t')

            for row in tsv:
                # Fetch ID, name and structure (SMILES)
                did = row[0]
                name = row[1]
                smiles = row[6]
                # Insert into SQL table
                self.sql.execute("INSERT INTO DRUGS VALUES(?,?,?)", (did, name, smiles))

        # Commit changes
        self.conn.commit()

    def get_drug_matches(self):
        '''
        Find variants that are connected to our drugs of interest
        And store them in the 'Variant' table
        '''

        print "Getting variants connected to drugs... (- w-)b"
        # Remove and re-create table base (initiate empty table)
        self.remake_table("drugvars")
        # Get all drugs from the drugs table
        self.sql.execute('''
                         SELECT DISTINCT DrugID
                         FROM Drugs
                         WHERE DrugID NOT IN (
                         SELECT DrugID FROM    DrugVars
                         )
                         ''')

        for (did,) in tqdm(self.sql.fetchall()):
            # Build URI to query API
            uri = \
                'https://api.pharmgkb.org/v1/report/connectedObjects/{}/Variant'.format(did)
            # Fetch results
            data = get_json(uri, self.authobj)
            if data is not None:
                print "data = ",
                print data
            # Fill INSERT SQL template with results
            sql = self.insert_sql("drugvars").render(json=data,
                                                     did=did)
            print "sql = " + sql
            # Execute insert statement
            time.sleep(2)
            self.sql.executescript(sql)
            # Commit changes
            self.conn.commit()

    def get_drug_vars(self):
        '''
        Get more information on the variant IDs we found in get_drug_matches
        and store in new table
        '''
        # Remove and re-create table base (initiate empty table)
        self.remake_table("variants")
        # Get template for filling 'Variant' table
        template = self.insert_sql("variants")
        # Get variant IDs found earlier
        self.sql.execute('''
                    SELECT DISTINCT VarID from DrugVars
                    WHERE VarID NOT IN (
                    SELECT DISTINCT VarID from Variants
                    );'''
                         )

        print "Getting more info on these variants..."
        # Store URI that will be queried
        uris = []

        # Loop through variant IDs
        for (varid,) in tqdm(self.sql.fetchall()):
            # Build URI for API querying
            uri = \
                'https://api.pharmgkb.org/v1/data/variant/{}?&view=max'.format(varid)
            uris.append(uri)
            # Query API and get results
            data = get_json(uri, self.authobj)
            # Fill template with resulting JSON file
            sql = template.render(json=data)
            # Execute SQL statement
            self.sql.executescript(sql)
            # Commit changes
            self.conn.commit()

    def get_gene_data(self):
        '''
        Fetches data on given gene IDs and
        inserts this information in new table
        :return:
        '''
        # Remove and re-create table base (initiate empty table)
        self.remake_table("genes")
        print 'Getting gene data... (/o*)'

        # Get all unique gene ids from the variant table
        self.sql.execute('''SELECT DISTINCT GeneID FROM Variants
        WHERE GeneID NOT IN (
            SELECT GeneID FROM Genes)''')

        # Loop through gene IDs
        for (gid,) in tqdm(self.sql.fetchall()):
            # Build URI for API querying
            uri = 'https://api.pharmgkb.org/v1/data/gene/{}?view=max'.format(gid)
            # Call API and fetch response (doesn't require authentication)
            data = urllib2.urlopen(uri)
            response = json.load(data)
            # Insert in 'Genes' table template
            try:
                sql = self.insert_sql("genes").render(json=response)
            except:
                pp.pprint(response)
                continue

            # Execute SQL statement
            self.sql.executescript(sql)
            # Commit changes
            self.conn.commit()
            # Sleep so server isn't overburdened
            time.sleep(3)

    def get_haplotypes(self):  # There is also a function with the same name defined in "patient.py"
        '''
        Finds known haplotypes for known Gene IDs
        and inserts in new table
        '''
        # Remove and re-create table base (initiate empty table)
        self.remake_table("haplotypes")
        # Define template for inserting in haplotype table
        template = self.insert_sql("haplotypes")
        # Get all unique gene IDs from the variant table
        self.sql.execute('''
        SELECT DISTINCT GeneID FROM Variants
        WHERE GeneID NOT IN (
        SELECT GeneID FROM Haplotypes)''')

        genes = self.sql.fetchall()

        # Go through results and fetch info connected to gene ID
        for (gid,) in tqdm(genes):
            # Build URI for API querying
            uri = 'https://api.pharmgkb.org/v1/data/haplotype?gene.accessionId={}&view=max'.format(
                gid)
            # Call API and fetch response
            data = get_json(uri, self.authobj)
            # Check if data was returned, if not go to next loop
            if not data:
                continue
            for response in data:
                # Insert in table template
                sql = template.render(json=response)
                # Execute SQL statement
                self.sql.executescript(sql)
            # Commit changes
            self.conn.commit()

    def get_hap_vars(self):
        """
        This function goes through variants connected to Haplotypes
        gets more information on them and adds their info to the 'Variant' table.
        """
        print 'Getting variants connected to haplotypes... ~(^_^)~'
        # Get RSIDs that are not present in the Variant table.
        self.sql.execute('''
                    SELECT DISTINCT VarName from HapVars
                    WHERE VarName LIKE "rs%" AND VarName NOT IN
                    (SELECT DISTINCT RSID FROM Variants)
                    ''')

        # Define template for inserting in variant table
        template = self.insert_sql("variants")

        # Loop through rsids
        for (rsid,) in tqdm(self.sql.fetchall()):
            # Build URI for API querying
            uri = \
                'https://api.pharmgkb.org/v1/data/variant/?symbol={}&view=max'.format(rsid)
            # Call API and fetch response
            data = get_json(uri, self.authobj)
            # Check if data was returned, if not go to next loop
            if data is None:
                continue
            # Insert in table template
            sql = template.render(json=data[0])
            # Execute SQL statement
            self.sql.executescript(sql)
            # Commit changes
            self.conn.commit()

    def get_non_rs(self):
        '''
        This function takes variants that DO NOT have a RSID
        and attempts to get the necessary information out of them.
        '''
        # Get non-RSID variants that are NOT HG38
        # (please adjust if you use hg38)

        self.sql.execute('''
                SELECT DISTINCT v.VarName, v.AltAllele, h.GeneID from HapVars v
                JOIN Haplotypes h ON v.HapID = h.HapID
                JOIN Genes g on h.GeneID = g.GeneID
                WHERE v.VarName LIKE "%chr%"
                OR v.VarName LIKE "chr%"
                AND v.VarName NOT LIKE "hg38"
                ORDER BY h.GeneID;
                ''')

        print 'Parsing non-rs variants...'

        # Define template for inserting in variant table
        template = self.insert_sql("othervars")

        # Loop through sql returned variant names
        for (rsid, alt, gid) in tqdm(self.sql.fetchall()):
            # Convert variant name to JSON format
            d = hg19conv(rsid, alt, gid)
            # Insert in table template
            sql = template.render(json=d)
            # Execute SQL statement
            self.sql.executescript(sql)
            # Commit changes
            self.conn.commit()

    def convert_indels(self):
        '''
        This function adjusts indel naming slightly to match GATK
        results. Example:
        deletion of *C* in *ACG*
        deletion in GATK:     REF: AC | ALT: A
        deletion in PharmGKB: REF:  C | ALT: del
        This is necessary otherwise no annotations return from
        the API calls later for these variants.
        '''

        # Define template for inserting in variant table
        template = self.insert_sql("indels")
        # Remove and re-create table base (initiate empty table)
        self.remake_table("indels")
        # Get necessary information for conversion from
        # the variant and variant location table
        self.sql.execute('''
                        SELECT DISTINCT
                        l.VarID, RefGenome, Chromosome, Start, End, RefAllele, AltPGKB
                        FROM LocPGKB l
                        JOIN Variants v
                        ON l.VarID = v.VarID
                        JOIN AltAlleles a
                        ON v.VarID = a.VarID
                        ''')

        print "Converting indels.. \(>w <)/"

        # Loop through results
        for (varid, genome, loc, start, end, ref, alt) in tqdm(self.sql.fetchall()):
            # Create json for template usage
            shifted = {"varid": varid,
                       "chromosome": loc,
                       "genome": genome,
                       "ref": ref,
                       "alt": alt,
                       "start": start,
                       "end": end}
            # Insertion or deletion?
            if alt == ref:
                continue

            if ref == "-":
                # Left shift position by 1
                shifted['start'] = start - 1
                shifted['end'] = end
                # Get reference nucleotide at that position
                prevbase = get_ref(loc, shifted['start'], shifted['start'])
                # Insertion scenario
                shifted['ref'] = prevbase
                shifted['alt'] = prevbase + alt

            elif alt == "-":
                # Left shift position by 1
                shifted['start'] = start - 1
                shifted['end'] = end
                # Get reference nucleotide at that position
                prevbase = get_ref(loc, shifted['start'], shifted['start'])
                # Deletion scenario A | -
                shifted['ref'] = prevbase + ref
                salt = prevbase
                shifted['alt'] = salt

            else:
                pass

            # Insert resulting JSON in template
            sql = template.render(alt=alt, json=shifted)
            # Execute sql query
            self.sql.executescript(sql)
            # Commit changes
            self.conn.commit()
# ------------------------------------------------------------------------------------------------------------

    def get_annotations(self):
        '''
        This function grabs annotations connected to drugs and variants in our drugs table.
        And stores their information in the Annotations table.
        '''
        # Remove and re-create table base (initiate empty table)
        self.remake_table("annotations")
        # Get known drug IDs from our table
        self.sql.execute('SELECT DISTINCT DrugID FROM DrugVars \
                                    WHERE DrugID NOT IN (SELECT DISTINCT DrugID FROM Annotations)')

        # Loop through drug IDs
        for (DrugID,) in tqdm(self.sql.fetchall()):
            # Build URI for API querying
            uri = \
                'https://api.pharmgkb.org/v1/data/clinicalAnnotation?relatedChemicals.accessionId={}' \
                .format(DrugID)
            # Call API and fetch response
            data = get_json(uri, self.authobj)
            # Insert JSON in template
            sql = self.insert_sql("annotations").render(json=data, DrugID=DrugID)
            # Execute insert statement
            self.sql.executescript(sql)
            # Commit changes
            self.conn.commit()

    def get_guidelines(self):
        '''
        Searches for drug dose guidelines connected to known drugs.
        Inserts results in Guidelines table.
        '''
        # Define template for later use
        template = self.insert_sql("guidelines")
        # Remove and re-create table base (initiate empty table)
        self.remake_table("guidelines")
        # Get known drug IDs.
        self.sql.execute('SELECT DISTINCT DrugID FROM DrugVars')

        # Loop through drug IDs.
        for (DrugID,) in tqdm(self.sql.fetchall()):
            # Build URI for API querying
            uri = 'https://api.pharmgkb.org/v1/data/guideline?relatedChemicals.accessionId={}'\
                .format(DrugID,)
            # Call API and fetch response
            data = get_json(uri, self.authobj)
            # Check if data was returned, if not go to next loop
            if data is None:
                continue
            # Insert JSON in template
            sql = template.render(json=data, did=DrugID)
            # Execute insert statement
            self.sql.executescript(sql)
        # Commit changes
        self.conn.commit()

    def get_guide_options(self):
        '''
        Gets Haplotype options for guidelines (not always present)
        CURRENTLY NON FUNCTIONAL
        '''
        # Define template for later use
        template = self.insert_sql("guideoptions")
        # Remove and re-create table base (initiate empty table)
        self.remake_table("guideoptions")
        # Get known guideline IDs.
        self.sql.execute("SELECT DISTINCT GuID, markdown from Guidelines;")

        # Loop through guideline IDs.
        for (guid, markdown) in tqdm(self.sql.fetchall()):
            # Build URI for API querying
            uri = "https://api.pharmgkb.org/v1/report/guideline/{}/options" \
                .format(guid)
            # Call API and fetch response
            data = get_json(uri, self.authobj)

            # Insert JSON in template
            # sql = template.render(guid = guid, json = data)

            # Execute insert statement
            # self.sql.executescript(sql)
        # Commit changes
        self.conn.commit()

    def bed_file(self):
        '''
        Creates bed file for subsetting .BAM files.
        This file has columns 'gene name', 'chromosome (CAREFUL w/ UCSC LABELING!!!)
        starting position and end position. To be used for running the GATK pipeline
        to get an interval GVCF file (so you don't need to create the whole GVCF).
        '''
        # Get locations of all genes involved in the database
        # ORDERING IS IMPORTANT!!
        self.sql.execute(
            'SELECT chromosome, start, end, genesymbol FROM genes ORDER BY chromosome ASC, start ASC')
        # lists for storing INT chromosomes (non-XY and mitocrondrial) and the rest
        linesINT = []
        linesETC = []

        print("here...")

        # Loop through genes
        for (chrom, start, end, name) in self.sql.fetchall():
            chrom = chrom.lstrip("chr")
            try:
                int(chrom)
                linesINT.append('\t'.join(map(str, (chrom, start, end, name))) + '\n')
            except:
                linesETC.append('\t'.join(map(str, (chrom, start, end, name))) + '\n')

        # Order by chromosome number
        linesINT.sort(key=lambda x: int(x.split("\t")[0]))
        # Join the two lists together
        linesALL = linesINT + linesETC
        # !!!!!!!!!!!UCSC SPECIFIC: DELETE IF NOT UCSC!!!!!!!!!!
        linesALLchr = ["chr" + line for line in linesALL]

        # Write to file :-)
        with open('PharmacogenomicGenes_PGKB_full.bed', 'w') as bed:
            bed.writelines(linesALLchr)

# -----------------------------------------------------------------
