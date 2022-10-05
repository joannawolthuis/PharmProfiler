#!/usr/bin/python
# -*- coding: utf-8 -*-

import requests
from requests_oauthlib import OAuth2Session
import getpass
import urllib
import urllib2
import ast
import re
import time

# ---------------------------------------------------------------------


def authenticate():  # This function has the same name as in "data.py"
    """
    This function creates an authenticating object for use in pharmgkb.
    Necessary for accessing clinical annotations and recommendations!
    :return:
    """
    print 'Authenticating...'

    email = raw_input('PharmGKB email adres: ')
    password = getpass.getpass()
    # email = raw_input('PGKB e-mail address: ')
    # password = getpass.getpass()
    req = {'email': email, 'password': password}

    # uri for authenticating with this file
    url = 'https://api.pharmgkb.org/v1/auth/oauthSignIn'
    # encode request including user info
    data = urllib.urlencode(req)
    # send request to specified url
    req = urllib2.Request(url, data)
    # read response
    response = urllib2.urlopen(req)
    token = response.read()
    # convert token to something that can be used for authentication
    token = ast.literal_eval(token)['accessToken']
    # create an authorized session and
    # return this session
    client = OAuth2Session(token=token, auto_refresh_url=url)
    client.access_token = token
    return client


def get_json(uri, client):
    """
    Uses uri and given authenticated client to get json data.
    :param uri: filled in uri refining query
    :param client: result of authenticate()
    :return:
    """
    resp = requests.head(uri)
    status = resp.status_code

    print resp
    print status
    print uri

    if status == 200:
        r = client.get(uri)
        print r
        data = r.json()
        return data

    else:
        return None


def seq_maker(rsidorder, reference, rsids):
    seq = ''

    for rsid in rsidorder:
        try:
            base = rsids[rsid]
        except KeyError:
            base = reference[rsid]
        seq += base
    return seq


def get_ref(loc, start, end):
    server = 'http://grch37.rest.ensembl.org'
    ext = '/sequence/region/human/{}:{}..{}?'.format(loc.lstrip('chr'),
                                                     start, end)
    r = requests.get(server + ext, headers={'Content-Type': 'text/plain'
                                            })

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    time.sleep(0.3)
    return r.text


def hg19conv(rsid, alt, gid):
    d = {}
    d['id'] = rsid
    d['gid'] = gid
    # find chromosome number
    d['loc'] = rsid.split(':')[0]
    # find location
    d['begin'] = int(rsid[rsid.find(':') + 1:rsid.find('(')])
    d['end'] = d['begin']
    # get reference position
    d['ref'] = get_ref(d['loc'], d['begin'], d['begin'])
    d['alt'] = alt

    if 'delGENE' in alt:
        return None

    if d['ref'] == alt:
        return None

    if len(d['ref']) == len(alt):
        d['muttype'] = 'snp'

    if 'ins' in alt:
        d['begin'] += 1
        d['ref'] = "-"
        d['alt'] = alt.replace('ins', "")
        d['muttype'] = 'in-del'
        d['end'] = d['begin'] + len(alt)

    elif 'del' in alt:
        d['alt'] = "-"
        d['muttype'] = 'in-del'
        d['end'] = d['begin'] + len(d['ref'])

    return d


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def deconstruct_guideline(markdown):
    markdown = markdown.split("\n")
    tex = ""
    # ------------
    nontable = []
    tables = {}
    cols = "XXXX"
    repline = "reference allele at all positions"
    for line in markdown:
        if repline in line:
            line = line.replace(repline, "Reference | Reference")
        if "|" in line:
            spl_line = [item for item in line.split("|") if item.strip() != ""]
            cols = len(spl_line) * "X"
            # ----------------------------------
            try:
                if "---" in spl_line[0]:
                    tables[cols].append("\hline")
                else:
                    tables[cols].append("&".join(spl_line))
            except:
                tables[cols] = ["&".join(["\\textbf{%s}" % item for item in spl_line])]
        # -------------------------------------
        else:
            nontable.append(line)

    return {'nontable': nontable, "tables": tables}
