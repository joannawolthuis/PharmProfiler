#!/usr/bin/python
# -*- coding: utf-8 -*-

from modules.pgkb_functions import *
from tqdm import tqdm
import pprint as pp
from db import Database
import os

# --------------------------------------------------------------------------


class Interpreter(Database):

    def __init__(self, patObj):

        self.sql = patObj.sql
        self.conn = patObj.conn
        self.tempfolder = patObj.tempfolder
        self.templateEnv = patObj.templateEnv
        self.advice = {}
        self.authobj = authenticate()
        self.conn.commit()

    def get_genotype(self):

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
                pat.RefAllele as PatRef,
                pat.AltAllele as PatAlt,
                CallNum,
                locv.start as start

                FROM hapvars h
                   JOIN
                   variants v ON v.rsid = h.VarName
                   JOIN
                   locvcf locv ON locv.varid = v.varid
                   JOIN
                   pat.Variants pat ON pat.start = locv.start
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
        self.sql.executescript('''
                        DROP TABLE IF EXISTS pat.Genotypes;
                        CREATE TABLE IF NOT EXISTS pat.Genotypes
                        (GeneID text,
                        New1_1 text, New1_2 text, New1_3 text,
                        New2_1 text, New2_2 text, New2_3 text,
                        Old1_1 text, Old1_2 text, Old1_3 text,
                        Old2_1 text, Old2_2 text, Old2_3 text);
        ''')

        self.sql.execute('''
        SELECT DISTINCT h.GeneID, g.GeneSymbol from Haplotypes h
        JOIN Genes g on g.GeneID = h.GeneID''')

        genotype = {}

        for (gid, genesymbol) in self.sql.fetchall():

            print genesymbol
            self.sql.executescript('''
            DROP VIEW IF EXISTS CurView;
            create temp view CurView AS
            SELECT DISTINCT
            h.GeneID as GeneID,
            p.HapID,
            p.Distance1 as new1, p.Distance2 as new2,
            po.score1 as old1, po.score2 as old2,
            po.MatchLen,
            h.hgvs, starname
            FROM pat.Haplotypes p
            JOIN pat.Haplotypes_OLD po
            on p.hapid = po.hapid
            JOIN Haplotypes h
            on h.HapID = P.HapID
            ''')

            haplotypes = self.sql.fetchall()

            self.sql.execute('''
    SELECT DISTINCT starname, new1, new2, old1, old2 FROM CurView
    WHERE GeneID = ?
    and hgvs LIKE "%=%"
                ''', (gid,))

            refhap = {}

            for (starname, new1, new2, old1, old2) in self.sql.fetchall():
                refhap['starname'] = starname
                refhap['new1'] = float(new1)
                refhap['new2'] = float(new2)
                refhap['old1'] = 0.0
                refhap['old2'] = 0.0

            allelesNEW = {1: [], 2: []}
            allelesOLD = {1: [], 2: []}
            allelesNEW = {1: [], 2: []}
            allelesOLD = {1: [], 2: []}

            # NEW STYLE
            for i in range(1, 3):
                self.sql.execute('''
                SELECT DISTINCT starname, new{} FROM CurView
                WHERE GeneID = ?
                ORDER BY new{} ASC
                '''.format(i, i), (gid,))

                # shuffle for ref at the top
                for (starname, score) in self.sql.fetchall():
                    if float(score) == refhap["new{}".format(i)] and starname != refhap['starname'] and starname not in allelesNEW[i]:
                        allelesNEW[i].append(refhap['starname'])
                        allelesNEW[i].append(starname)
                    else:
                        allelesNEW[i].append(starname)

                # OLD STYLE
                self.sql.execute('''
                SELECT DISTINCT starname, old{}, matchLen FROM CurView
                WHERE GeneID = ?
                ORDER BY old{} DESC, MatchLen DESC
                '''.format(i, i), (gid,))

                for (starname, score, mlen) in self.sql.fetchall():
                    if float(score) == 0.0 and starname != refhap['starname'] and starname not in allelesOLD[i]:
                        allelesOLD[i].append(refhap['starname'])
                        allelesOLD[i].append(starname)
                    else:
                        allelesOLD[i].append(starname)

            if len(allelesNEW[1]) == 0:
                continue
            elif len(allelesNEW[1]) >= 3:
                item = (gid, allelesNEW[1][0], allelesNEW[1][1], allelesNEW[1][2],
                        allelesNEW[2][0], allelesNEW[2][1], allelesNEW[2][2],
                        allelesOLD[1][0], allelesOLD[1][1], allelesOLD[1][2],
                        allelesOLD[2][0], allelesOLD[2][1], allelesOLD[2][2])
            elif len(allelesNEW[1]) == 2:
                item = (gid, allelesNEW[1][0], allelesNEW[1][1], "n/a",
                        allelesNEW[2][0], allelesNEW[2][1], "n/a",
                        allelesOLD[1][0], allelesOLD[1][1], "n/a",
                        allelesOLD[2][0], allelesOLD[2][1], "n/a")
            elif len(allelesNEW[1]) == 1:
                item = (gid, allelesNEW[1][0], "n/a", "n/a",
                        allelesNEW[2][0], "n/a", "n/a",
                        allelesOLD[1][0], "n/a", "n/a",
                        allelesOLD[2][0], "n/a", "n/a")
            else:
                continue
            print item
            self.sql.execute("INSERT INTO pat.Genotypes VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)", item)
            self.conn.commit()

    def find_guidelines(self):

        self.remake_table("patguidelines")
        # For each drug:
        self.sql.execute("SELECT DISTINCT GuID from Guidelines")

        for (guid,) in tqdm(self.sql.fetchall()):
            print guid
            # fetch involved genes
            self.sql.execute(
                "SELECT DISTINCT o.GeneID, g.GeneSymbol FROM Guidelines o JOIN Genes g on g.GeneID = o.GeneID WHERE GuID = ?", (guid,))
            genes = self.sql.fetchall()
            guideGenes = {}
            print genes

            if len(genes) == 0:
                continue
            # Fetch involved haplotypes and scores

            genotype = {}

            for (gid, genesymbol) in genes:
                self.sql.execute(
                    "SELECT DISTINCT Starname FROM GuideOptions o JOIN Genes g on g.GeneSymbol = o.GeneSymbol WHERE GuID = ?", (guid,))

                options = [starname[0] for starname in self.sql.fetchall()]
                print options
                self.sql.execute('''
                SELECT DISTINCT
                New1_1, New2_1,
                Old1_1, Old2_1
                FROM pat.Genotypes
                WHERE GeneID = ?
                ''', (genesymbol,))

                scores = self.sql.fetchone()
                print scores

                if scores is None:
                    continue

                guideGenes[gid] = {}
                guideGenes[gid]['new'] = "{}:{}/{}".format(genesymbol, scores[0], scores[1])
                guideGenes[gid]['old'] = "{}:{}/{}".format(genesymbol, scores[2], scores[3])

        # ----------------------------------------------------------------------
                string_genotype_NEW = ";".join([guideGenes[gid]['new']
                                                for gid in guideGenes.keys()])
                string_genotype_OLD = ";".join([guideGenes[gid]['old']
                                                for gid in guideGenes.keys()])

                print "NEW:", string_genotype_NEW
                print "OLD:", string_genotype_OLD

                # Find matching advice
                string_genotype = string_genotype_OLD
                uri = "https://api.pharmgkb.org/v1/report/guideline/{}/annotations?genotypes={}".format(
                    guid, string_genotype)
                data = get_json(uri, self.authobj)

                if data is None:
                    uri = "https://api.pharmgkb.org/v1/report/guideline/{}/annotations?genotypes={}".format(
                        guid, string_genotype.replace("1A", "1"))
                    data = get_json(uri, self.authobj)

                sql = self.insert_sql("patguidelines").render(
                    guid=guid, genotype=string_genotype, json=data)
                print sql
                self.sql.executescript(sql)
                self.conn.commit()
                # Save to advice table 'PatGuidelines' (DrugID, GeneID, Category(Metabolizer type), Advice)

    def make_annotation(self):

        self.remake_table("patannotations")
        # only when there is no haplotype available?
        self.sql.execute('''
                            SELECT DISTINCT
                                    a.AnID,
                                    VarName,
                                    VarID,
                                    MutType,
                                    RefPGKB,
                                    AltPGKB,
                                    RefVCF,
                                    AltVCF,
                                    PatRef,
                                    PatAlt,
                                    CallNum,
                                    Start
                            FROM Overview
                            JOIN Annotations a
                            on a.VarHapID = VarID
                            WHERE RefVCF = PatRef
                            AND a.AnID NOT IN
                            (SELECT DISTINCT AnID FROM pat.Annotations)''')

        print "Annotating SNPs... /(* ` ^ `*/)"

        for (anid, rsid, varid, muttype, refP, altP, refV, altV, ref, alt, call, start) in tqdm(self.sql.fetchall()):
            uri = 'https://api.pharmgkb.org/v1/data/clinicalAnnotation/{}?view=max'.format(anid)
            data = get_json(uri, self.authobj)
            # --------------------------------------
            if muttype == "snp":
                call = call.replace("/", "")
                allele = call.replace("0", ref)
                allele = allele.replace("1", alt)
                revAllele = allele[::-1]

            elif muttype != "snp":
                allele = call.replace("0", refP.replace("-", "del"))
                allele = allele.replace("1", altP.replace("-", "del"))
                s_allele = allele.split("/")
                revAllele = s_allele[1] + "/" + s_allele[0]
            # --------------------------------------
            sql = self.insert_sql("patannotations").render(
                json=data, revallele=revAllele, patallele=allele)
            print sql
            # try:
            self.sql.executescript(sql)
            # except:
            #    print sql
            #    continue
        self.conn.commit()

# ---------------------------- NEXT STEP: ReportMaker --------------------------------
