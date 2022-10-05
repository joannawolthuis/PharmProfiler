#!/usr/bin/python
# -*- coding: utf-8 -*-

from db import *
from data import *
import os
import vcf
from tqdm import tqdm
from modules.pgkb_functions import *
import subprocess as s
import pprint as pp
from Bio import Phylo, AlignIO
from collections import OrderedDict
from interpreter import Interpreter

# ---------------------------------------------------------------------


class Patient(Database):
    '''
    This class imports a design VCF and an input VCF.
    '''

    def __init__(self, dbname, f):
        """
        f = given VCF input file
        GetDesign imports the test design and matches rs values to
        positions
        Import imports the input vcf
        GetIDs finds changed positions and matches them to the design.
        """

        print 'Initiating patient...'

        self.path = os.path.dirname(__file__)
        dbfolder = os.path.join(self.path, 'db')
        dbpath = os.path.join(dbfolder, '%s.db' % dbname)
        Database.__init__(self, dbpath)
        self.fname = os.path.basename(f).rstrip(".g.vcf.gz")
        self.patDb = f.replace(".g.vcf.gz", ".db")

        print self.patDb

        self.sql.execute("attach '{}' as pat".format(self.patDb))
# -------------------------------------------------------------
        self.f = f
        self.reader = vcf.Reader(open(f, 'r'))
# --------------------------------- SNP parsing ----------------------------------------

    def get_positions(self):
        """
        Gets patient variables, reads into patientvars table
        :return:
        """
        self.remake_table("patvars")
        # create list of positions to localize in patient
        self.sql.execute('''
            SELECT DISTINCT l.Chromosome, l.Start, l.End, l.RefAllele, v.MutType, v.VarID
            FROM LocVCF l
            JOIN Variants v
            ON v.VarID = l.VarID
            ''')

        positions = self.sql.fetchall()

        for (loc, start, end, ref, muttype, varid) in tqdm(positions, desc="Fetching patient variants o(* ^ *o)"):
            loc = "chr" + str(loc.lstrip("chr"))
            print loc
            try:
                records = self.reader.fetch(loc, start - 1, end=end)
            except:
                print "Couldn't find {} in gvcf".format(varid)
                continue

            for record in records:  # doctest: +SKIP
                sql = self.insert_sql("patvars").render(record=record)
                self.sql.executescript(sql)

# ---------------------------- Haplotype parsing -------------------------------

    def get_haplotypes(self):  # There is also a function with the same name defined in "data.py"
        """
        Matches rsids per gene to known haplotypes for that gene.
        Results are stored in a list...
        :return:
        """
        self.remake_table('pathaplotypes')
        self.conn.commit()
        self.sql.execute("SELECT DISTINCT GeneID from Haplotypes")
        gids = [tup[0] for tup in self.sql.fetchall()]
        self.sql.executescript('''DROP TABLE IF EXISTS pat.Haplotypes_OLD;
                                                CREATE TABLE pat.Haplotypes_OLD(hapid text, score1 text, score2 text, matchLen text);''')

        # create view with everything
        self.sql.executescript('''
                    DROP VIEW IF EXISTS Overview;

                    CREATE TEMP VIEW Overview AS
                    SELECT DISTINCT

                    hap.GeneID as GeneID,
                    hap.HapID as HapID,
                    hap.HGVS as HGVS,
                    hap.starname as starname,
                    VarName, MutType,
                    v.VarID as VarID,
                    h.AltAllele as HapAllele,
                    locb.RefAllele as RefPGKB,
                    a.AltPGKB as AltPGKB,
                    locv.RefAllele as RefVCF,
                    a.AltVCF as AltVCF,
                    p.RefAllele as PatRef,
                    p.AltAllele as PatAlt,
                    CallNum,
                    locv.start as start

                    FROM hapvars h
                       JOIN
                       variants v ON v.rsid = h.VarName
                       JOIN
                       locvcf locv ON locv.varid = v.varid
                       JOIN
                       pat.Variants p ON p.start = locv.start
                       JOIN
                       haplotypes hap ON hap.hapid = h.HapID
                       JOIN
                       genes gen ON gen.geneid = hap.geneid
                       JOIN
                       locpgkb locb ON locb.varid = locv.varid
                       JOIN
                       altalleles a on a.VarID = locv.VarID
                       ORDER BY locv.start asc;
                            ''')

        # get list of all gids
        print "Haplotyping patient... (\\' n')\\*scribble*"
        scoreOverview = {}

        for gid in tqdm(gids):
            varValues = OrderedDict([])

            # First: collect SNPs and patient positions at those SNPs
            self.sql.execute('''
                    SELECT DISTINCT hapid, VarName, HapAllele, PatAlt, RefPGKB, AltPGKB, RefVCF, AltVCF, CallNum, MutType
                    FROM Overview
                    WHERE GeneID = ?
                    AND HGVS LIKE "%[=]%"
                    AND HapAllele NOT LIKE "%(%"
                    AND HapAllele NOT LIKE "%[%"
                    ORDER BY Start ASC
                    ''', (gid,))

            refRsids = self.sql.fetchall()

            rsidorder = list(set([rsid
                                  for (hapid, rsid, hapallele, patAlt, refp, altp, refv, altv, CallNum, MutType)
                                  in refRsids]))

            if len(rsidorder) == 0:
                continue

            snps = {}
            indels = {}
            patrsids_hom = {}
            patrsids_het = {}

            for (hapid, rsid, hapallele, patAlt, refp, altp, refv, altv, CallNum, MutType) in refRsids:
                refid = hapid
                # SNPs
                if "[" not in hapallele and "(" not in hapallele:
                    if MutType == "snp":
                        snps[rsid] = hapallele

                    elif MutType != "snp":
                        if hapallele == refp or hapallele[::-1] == refp:
                            indels[rsid] = refv
                        elif hapallele == altp or hapallele[::-1] == altp:
                            indels[rsid] = altv
                        else:
                            print hapid, rsid, hapallele, MutType, CallNum, refp, altp, refv, altv
                            continue

                    if CallNum == "1/1":
                        patrsids_hom[rsid] = patAlt
                        patrsids_het[rsid] = patAlt
                    if CallNum == "0/1":
                        patrsids_het[rsid] = patAlt

                else:
                    print "cannot parse", rsid, hapallele, MutType
                    continue
                # GET PATIENT VARS

            refid = refRsids[0][0]
            scoreOverview[refid] = {"al1": 0.0, "al2": 0.0, "hapLen": 0}

            # Join these three dictionaries and feed to seq_maker
            varValues[refid] = merge_dicts(snps, indels)
            varValues['a1'] = merge_dicts(varValues[refid], patrsids_hom)
            varValues['a2'] = merge_dicts(varValues[refid], patrsids_het)

            # get list of all hapids for this gene
            self.sql.execute('''
                    SELECT DISTINCT HapID, starname
                    from Overview
                    where GeneID = ?
                    and hgvs not like "%=%"
                    ''',
                             (gid,))

            selection = self.sql.fetchall()

            for (hapid, starname) in selection:
                print hapid, starname
                scoreOverview[hapid] = {}
                # get haplotype alleles and create complete dictionary
                self.sql.execute('''SELECT VarName, AltVCF
                            from Overview
                            where HapID = ?
                            ''', (hapid,))

                haprsids = {rsid: alt for (rsid, alt) in self.sql.fetchall()}
                self.sql.execute("INSERT INTO pat.Haplotypes_OLD VALUES(?,?,?,?)", (refid, 0, 0, 0))

                if len(haprsids.items()) == 0:
                    print "no rsids found for", hapid
                    continue

                else:
                    varValues[hapid] = dict(varValues[refid], **haprsids)
                    uniques_dct = dict(set(varValues[refid].items()) -
                                       set(varValues[hapid].items()))

                    if len(uniques_dct.items()) == 0:
                        print "no unique rsids found for", hapid
                        continue

                    uniques_al1 = dict(set(varValues[refid].items()) - set(varValues['a1'].items()))
                    uniques_al2 = dict(set(varValues[refid].items()) - set(varValues['a2'].items()))

                    # calculate match scores old way
                    scores = {}
                    hapLen = len(haprsids.keys())

                    shared_al1 = dict(set(uniques_dct.items()) & set(uniques_al1.items()))
                    shared_al2 = dict(set(uniques_dct.items()) & set(uniques_al2.items()))

                    match_score1 = float(len(shared_al1.keys())) / float(len(uniques_dct.keys()))
                    match_score2 = float(len(shared_al2.keys())) / float(len(uniques_dct.keys()))

                    self.sql.execute("INSERT INTO pat.Haplotypes_OLD VALUES(?,?,?,?)",
                                     (hapid, match_score1, match_score2, hapLen))
                    self.conn.commit()

# -------------------------------------------------------------------------------------------------------------------

            output = "/output/alignments/"
            fn = self.path + output + self.fname + "_" + gid + "_aln.fasta"
            prev_seqs = []

            with open(fn, "w") as f:
                refseq = seq_maker(rsidorder, varValues[refid], varValues[refid])
                prev_seqs.append(refseq)
                f.write(">{}\n{}\n".format(refid, refseq))

                for var, values in varValues.items():
                    seq = seq_maker(rsidorder, varValues[refid], values)

                    if seq not in prev_seqs or var == "a1" or var == "a2":
                        f.write(">{}\n{}\n".format(var, seq))

                    if var != "a1" and var != "a2":
                        prev_seqs.append(seq)

    def hap_scorer(self):
        output = "/output/alignments/"
        self.sql.execute("SELECT DISTINCT GeneID FROM Haplotypes")

        for (gid,) in tqdm(self.sql.fetchall(), desc="Scoring haplotypes..."):
            self.sql.execute(
                "SELECT DISTINCT HapID FROM Haplotypes WHERE hgvs LIKE '%[=]%' AND GeneID = ?", (gid,))
            refids = [hapid for (hapid,) in self.sql.fetchall()]
            fn = self.path + output + self.fname + "_" + gid + "_aln.fasta"

            # phylogenetic tree
            of = fn.replace("alignments/", "alignments/aligned/")
            tn1 = of.strip(".fasta")+".dnd"
            tn2 = of.strip(".fasta")+"_rooted.dnd"

            try:
                with open(fn, "r") as f:
                    s.check_call("{}/plugins/clustalo -i {} -o {} -t DNA --dealign --force"
                                 .format(self.path, fn, of),  shell=True)

            except:
                continue

            s.check_call(
                "{}/plugins/FastTree -quiet -nopr -gtr -nt {} > {}".format(self.path, of, tn1), shell=True)
            tree = Phylo.read(tn1, "newick")
            names = []

            for refid in refids:
                try:
                    tree.root_with_outgroup({refid})

                except:
                    continue

            # add tree to db
            Phylo.write(tree, tn2, "newick")

            for clade in tree.find_clades():
                if clade.name and clade.name not in names:
                    names.append(clade.name)

            for hap in names:
                distances = {"hapid": hap}
                if "a1" in hap or "a2" in hap:
                    continue

                else:
                    for i in range(1, 3):
                        dist = tree.distance(target1="a%i" % i, target2=hap)
                        distances["al%i" % i] = dist

                    sql = self.insert_sql("pathaplotypes").render(json=distances, gid=gid)
                    print sql
                    self.sql.executescript(sql)
                    self.conn.commit()
        self.conn.commit()

    def Interpret(self):
        # Check for reference alleles
        i = Interpreter(self)
        # i.Genotyper()
        i.make_annotation()
