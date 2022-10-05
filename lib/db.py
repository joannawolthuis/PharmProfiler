#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3
import os
import re
import glob
from jinja2 import Template, FileSystemLoader, Environment
from modules.pgkb_functions import *

# --------------------------------------------------------------------------


class Database(object):
    '''
 Database class. Initializes database and takes care of doing big changes.
    '''

    def __init__(self, dbpath):

        path = os.path.dirname(dbpath)
        print path
        self.tempfolder = path.replace("db", "templates")
        self.conn = sqlite3.connect(dbpath)  # connect to db- if it doesn't exist, create it
        self.sql = self.conn.cursor()  # cursor for sqlite3, used to do things in database
        templateLoader = FileSystemLoader(searchpath=self.tempfolder)
        self.templateEnv = Environment(loader=templateLoader)

    def load_sql(self, path):
        sql = ''

        with open(path, 'r') as f:
            for line in f.readlines():
                sql += re.sub('[\n\t]', '', line)
        return sql

    def get_templates(self, template):
        TEMPLATE_FILE = template
        template = self.templateEnv.get_template(TEMPLATE_FILE)
        return template

    def insert_sql(self, tabname):
        TEMPLATE_FILE = tabname + '.ins'
        template = self.templateEnv.get_template(TEMPLATE_FILE)
        return template

    def set_defaults(self):  # I could not find any place where this function is used?
        dropTables = glob.glob(os.path.join(self.tempfolder, "*.rm"))
        createTables = glob.glob(os.path.join(self.tempfolder, "*.tab"))

        for template in dropTables + createTables:
            sql = self.load_sql(template)

        self.sql.executescript(sql)
        self.conn.commit()
        # test

    def remove_table(self, tabnames):
        if type(tabnames) != list:
            tabnames = [tabnames]
        for tabname in tabnames:
            sql = self.load_sql(os.path.join(self.tempfolder, tabname + '.rm'))
            print sql
            self.sql.executescript(sql)
            self.conn.commit()

    def create_table(self, tabnames):
        if type(tabnames) != list:
            tabnames = [tabnames]
        for tabname in tabnames:
            sql = self.load_sql(os.path.join(self.tempfolder, tabname + '.tab'))
            print sql
            self.sql.executescript(sql)
            self.conn.commit()

    def remake_table(self, tabnames):
        self.remove_table(tabnames)
        self.create_table(tabnames)

# -----------------------------------------------------
