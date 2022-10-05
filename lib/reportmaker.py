#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3
from jinja2 import Template, FileSystemLoader, Environment
import os
from tqdm import tqdm
import sys
import pprint as pp
from db import Database
from patient import Patient
from collections import OrderedDict
from modules.pgkb_functions import deconstruct_guideline
import re


def tex_escape(text):
    """
        :param text: a plain text message
        :return: the message escaped to appear correctly in LaTeX
    """
    conv = {
        '%': r'\%',
        '$': r'\$',
        '#': r'\#',
        '_': r'',
        '^': r'-',
        "^a^": "",
        "^b^": "",
        "^c^": "",
        "^d^": "",
        u"\u0001": ""
    }
    regex = re.compile('|'.join(re.escape(unicode(key))
                                for key in sorted(conv.keys(), key=lambda item: - len(item))))
    return regex.sub(lambda match: conv[match.group()], text)


def convert_latex(obj):
    conv = OrderedDict([
        ("%", "\%"),
        ("_", ""),
        (";", ","),
        ("[", "("),
        ("]", ")"),
        ("^a^", ""),
        ("^b^", ""),
        ("^c^", ""),
        ("^d^", ""),
        ("&gt;", ">"),
        ("&lt;", "<"),
        (u"\u0001", ""),
        ("&#34,", "\""),
        ("#", "\#"),
    ])
    for before, after in conv.items():
        obj = obj.replace(before, after)
    return obj


def change_dict_naming_convention(d, convert_function):
    """
    Recursivly goes through the dictionnary obj and replaces keys with the convert function.
    """
    new = {}
    if isinstance(d, unicode) or isinstance(d, str):
        new = convert_function(d)
        return new
    else:
        for k, v in d.items():
            new_v = v
            if isinstance(v, dict):
                new_v = change_dict_naming_convention(v, convert_function)
            elif isinstance(v, list):
                new_v = list()
                for x in v:
                    new_v.append(change_dict_naming_convention(x, convert_function))
            elif isinstance(v, unicode) or isinstance(v, str):
                new_v = convert_function(v)
            new[convert_function(k)] = new_v
    return new


# Needed for document:
JSON_TEMPLATE_DESIGN = '''
# --------------------------------------------------------------------------------------

{samplename, haplotable:[ {symbol:{ new1:{1, 2, 3}, new2{1, 2, 3}, old1{1, 2, 3}, old2{1, 2, 3}, tree)}, ...],
                        annotable:[ (drug, gene, guidelinelink, Ann1-2, Ann3-4), ...],
                        drugs:[drugname:{drugDesc, genes:[{symbol, hapNEW, hapOLD, guideline, annotations:[:{1A,
                            1B:{patAllele, rsid, text},
                            2A, 2B, 3, 4}, ...]
                        ]}
                            }
                        } , ...]
}
'''
# ============================================================


class ReportMaker(Database):

    def __init__(self, patObj, outfile):
        self.sn = patObj.fname.split("_")[0].replace("BLOOD", "").replace("blood", "")
        print self.sn
        self.sql = patObj.sql
        self.conn = patObj.conn
        self.tempfolder = patObj.tempfolder
        self.templateEnv = patObj.templateEnv
        print "Initiating Report..."
        self.path = os.path.dirname(__file__)
        self.outfile = outfile
    # ==========================================================
        # AnnotationView
        self.template = self.get_templates("latex/python_template.tex")

    def make_json(self):
        # gather data for tables
        JSON = {"samplename": self.sn,
                "haplotable": [],
                "annotable": [],
                "drugs": {}
                }
        self.sql.execute("select name from pat.sqlite_master where type = 'table';")
        print self.sql.fetchall()
        self.sql.execute('''
                    select distinct g.genesymbol, g.geneid,
                    new1_1, new1_2, new1_3,
                    new2_1, new2_2, new2_3,
                    old1_1, old1_2, old1_3,
                    old2_1, old2_2, old2_3
                    from
                    pat.genotypes p
                    join genes g
                    on g.GeneID = p.GeneID
                    order by g.genesymbol
                                ''')

        for (symbol, gid, n11, n12, n13, n21, n22, n23, o11, o12, o13, o21, o22, o23) in tqdm(self.sql.fetchall(), desc="Creating haplotype table..."):
            entry = {"symbol": symbol,
                     "new1": {"1": n11, "2": n12, "3": n13},
                     "new2": {"1": n21, "2": n22, "3": n23},
                     "old1": {"1": o11, "2": o12, "3": o13},
                     "old2": {"1": o21, "2": o22, "3": o23}}
            JSON['haplotable'].append(entry)

#                tree = ""
#                with open(self.path + "/output/alignments/aligned/" + gid + "_aln_rooted.dnd", "rb") as f:
#                    self.sql.execute('SELECT DISTINCT hapid, starname FROM Haplotypes WHERE GeneID = ?', (gid,))
#                    hapconvert = {hapid:starname for (hapid, starname) in self.sql.fetchall()}
#                    content = f.readlines()
#                    for line in content:
#                        for hapid, starname in hapconvert.items():
#                            starname = starname.replace(",", "")
#                            if "(" in starname:
#                                starname = starname.split(" (")[0]
#                            line = line.replace(hapid, starname).replace("a1", "{}_hap1".format(self.sn)).replace("a2", "{}_hap2".format(self.sn))
#                        tree += line.replace(" ", "").rstrip("\n")
#                    entry['tree'] = tree

            # ------------------------------

        # get gene-drug pairs
        self.sql.executescript('''
                                DROP VIEW IF EXISTS JsonView;
                                CREATE TEMP VIEW JsonView AS
                                select distinct
                                g.genesymbol as symbol,
                                g.geneid as geneid,
                                g.genename as genename,
                                d.chemname as drugname,
                                d.drugid as drugid,
                                l.guid as guid,
                                l.markdown as markdown,
                                l.summary as summary,
                                a.loe as loe,
                                v.rsid as rsid,
                                p.patallele as allele,
                                p.phenotype as phenotype
                                from genes g
                                join variants v
                                on v.geneid = g.geneid
                                left join annotations a
                                on a.varhapid = v.varid
                                join drugs d
                                on a.drugid = d.drugid
                                join pat.annotations p
                                on p.anid = a.anid
                                left join guidelines l
                                on l.geneid = g.geneid
                                and l.drugid = d.drugid
                                order by d.chemname asc;
                                    ''')

        self.sql.execute('''
        SELECT DISTINCT symbol, geneid, genename, drugname, drugid, guid from JsonView
        ''')
        for (symbol, gid, genename, drugname, did, guid) in tqdm(self.sql.fetchall(), desc="Generating report..."):
            if drugname not in JSON['drugs'].keys():
                JSON['drugs'][drugname] = {"genes": []}
                hasGuide = False
                hasAnno12 = False
                hasAnno34 = False

            # get annotation amounts
            self.sql.execute('''
                    select distinct LoE, count(*) from JsonView
                    where drugid = ?
                    group by LoE;
                    ''',  (did,))

            for (loe, amount) in self.sql.fetchall():
                cat1 = 0
                cat2 = 0
                if "1" in loe or "2" in loe:
                    cat1 += amount
                elif "2" in loe or "3" in loe:
                    cat2 += amount
                else:
                    pass

            # entry for summary table 2
            annoEntry = {"drugname": drugname,
                         "annCount":
                         {"1-2": cat1, "3-4": cat2}
                         }
            mainEntry = {}
            if annoEntry not in JSON['annotable']:
                JSON['annotable'].append(annoEntry)

            # ------ get more data on guidelines and annotations -----
            geneEntry = {
                "genesymbol": symbol,
                "genename": genename,
            }

            for hap in JSON['haplotable']:
                if hap["symbol"] != symbol:
                    continue
                else:
                    geneEntry['hapNEW'] = "{}/{}".format(hap['new1']["1"], hap['new2']["1"])
                    geneEntry['hapOLD'] = "{}/{}".format(hap['old1']["1"], hap['old2']["1"])

            # ^ this goes in json[drugs][DRUGNAME][genes] through append to list ^

            # --- get the annotation data ---
            self.sql.execute('''
                        select distinct
                        loe, allele,
                        rsid, phenotype,
                        guid, markdown, summary
                        from JsonView
                        where geneid = ?
                        and drugid = ?
                        ''', (gid, did))

            for (loe, allele, rsid, phenotype, guid, markdown, summary) in self.sql.fetchall():
                # get guideline info
                if guid != None:
                    hasGuide = True
                    geneEntry["guideline"] = {}
                    geneEntry['guideline']['guid'] = guid
                    geneEntry['guideline']['summary'] = summary
                    geneEntry['guideline']['text'] = deconstruct_guideline(markdown)['nontable']
                    geneEntry['guideline']['tables'] = deconstruct_guideline(markdown)['tables']

                # get annotation info
                if phenotype != None:
                    if "1" in loe or "2" in loe:
                        hasAnno12 = True
                    if "3" in loe or "4" in loe:
                        hasAnno34 = True
                    try:
                        g = geneEntry['annotations']
                    except:
                        geneEntry["annotations"] = {"1A": [], "1B": [],
                                                    "2A": [], "2B": [], "3": [], "4": []}
                    annEntry = {"rsid": rsid, "patAllele": allele, "phenotype": phenotype}
                if annEntry not in geneEntry['annotations'][loe]:
                    geneEntry['annotations'][loe].append(annEntry)
            if geneEntry not in JSON['drugs'][drugname]['genes']:
                JSON['drugs'][drugname]['genes'].append(geneEntry)
                JSON['drugs'][drugname]['hasGuide'] = hasGuide
                JSON['drugs'][drugname]['hasAnno12'] = hasAnno12
                JSON['drugs'][drugname]['hasAnno34'] = hasAnno34

        from HTMLParser import HTMLParser
        h = HTMLParser()

        self.JSON = change_dict_naming_convention(JSON, h.unescape)

# -----------------------------------------------------------------------------------------

    def make_report(self):
        reportText = self.template.render(sampleName=self.sn, json=self.JSON)
        with open(self.outfile + self.sn + ".tex", "wb") as f:
            f.write(tex_escape(reportText.encode("utf-8")))
