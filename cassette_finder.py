# Author: Arkadiy Garber
# !/bin/sh
from collections import defaultdict
import re
import os
import textwrap
import argparse
import urllib.request
import ssl
from urllib.error import HTTPError

gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1)

parser = argparse.ArgumentParser(
    prog="cassette_finder.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    Program for locating gene cassettes and clusters;
    Developed by Arkadiy Garber; University of Delaware, Geological Sciences
    Please send comments and inquiries to arkg@udel.edu
    '''))

parser.add_argument('genes', help='fasta file containing genes of interest (i.e. genes that are expected to be in '
                                  'gene cassette)')
parser.add_argument('db', help='Location of the the NCBI nr database')
parser.add_argument('-t', type=int, help='Number of threads to use (default=1)', default=1)
parser.add_argument('-e', type=float, help='E-value cut-off for gene homology (default=1E-5)', default="1E-5")
parser.add_argument('-g', type=int, help='Max gap for gene separation (i.e. maximum number of genes between genes '
                                         'of interest) (default=40)', default=40)
parser.add_argument('-n', type=int, help='Minimum number of genes you would like to see in cassette or cluster '
                                         '(default=2)', default=2)
args = parser.parse_args()


def DictSorter(dict):
    emptyList = []
    for i in dict.keys():
        emptyList.append(i)
    return sorted(emptyList)


def cluster(data, maxgap):
    data = sorted(data)
    # data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def fasta(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    count = 0
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[str(count) + " " + header] = seq
                count += 1
                header = i[1:]
                seq = ''
            else:
                header = i[1:]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


def fastaPlain(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                seq = ''
            else:
                header = i[1:]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


cwd = os.getcwd()


os.system("mkdir " + cwd + "/output")
os.system("blastp -query " + args.genes + " -db " +args.db + " -num_threads " + str(args.t) + " -out " +
          cwd + "/output/" + args.genes + ".blastResults -evalue " + str(args.e) + " -outfmt 6")


blast = open(cwd + "/output/" + args.genes + ".blastResults", "r")
bDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in blast:
    bList = i.split("\t")
    gene = bList[0]
    accessions = bList[1]
    gi = int(accessions.split("|")[1])
    bDict[gi] = bList


x = cluster(bDict.keys(), maxgap=40)
outfile = open(cwd + "/output/syntenies.csv", "w")
outfile.write("organism" + "," + "gene" + "," + "identification" + "," + "blast_e-value" + "\n")
for i in x:
    if len(i) >= int(args.n):
        print(i)
        try:
            fp = urllib.request.urlopen("https://www.ncbi.nlm.nih.gov/protein/%s" % i[0], context=gcontext)
            mybytes = fp.read()
            mystr = mybytes.decode("utf8")
            fp.close()
            lines = mystr.split('\n')
            for j in lines:
                if re.findall(r'<h1>', j):
                    j = j.split("[")
                    j = j[1].split("]")
                    org = j[0]
            for id in i:
                outfile.write(org + "," + bDict[id][0] + "," + bDict[id][1] + "," + bDict[id][10] + "\n")
        except HTTPError:
            pass
