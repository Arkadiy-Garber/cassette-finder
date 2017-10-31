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

parser.add_argument('-genes', help='fasta file containing genes of interest (i.e. genes that are expected to be in '
                                  'gene cassette)')
parser.add_argument('-db', help='Location of the the NCBI protein database or '
                               'folder containing translated genomes in fasta format')
parser.add_argument('-db_type', type=str, help='the type of database to search in (folder, nr, or refseq)', default="nr")
parser.add_argument('-output', type=str, help='name of output folder to write results in', default="output")
parser.add_argument('-t', type=int, help='Number of threads to use (default=1)', default=1)
parser.add_argument('-e', type=float, help='E-value cut-off for gene homology (default=1E-5)', default="1E-5")
parser.add_argument('-g', type=int, help='Max gap for gene separation (i.e. maximum number of genes between genes '
                                         'of interest) (default=40)', default=40)
parser.add_argument('-n', type=int, help='Minimum number of genes you would like to see in cassette or cluster '
                                         '(default=2)', default=2)
args = parser.parse_args()


def deStr(string):
    outNum = []
    for i in string:
        try:
            i = int(i)
            outNum.append(i)
        except ValueError:
            if i == ".":
                break
    stringNum = ''.join(map(str, outNum))
    Num = int(stringNum)
    return Num


def cluster2(data, maxgap):
    data = sorted(data)
    ls = []
    for i in data:
        ls.append(deStr(i))
    data = sorted(ls)
    # data.sort(key=int)
    try:
        groups = [[data[0]]]
        for x in data[1:]:
            if abs(x - groups[-1][-1]) <= maxgap:
                groups[-1].append(x)
            else:
                groups.append([x])
        return groups
    except IndexError:
        return "Unresolved program failure"


def cluster(data, maxgap):
    data = sorted(data)
    ls = []
    for i in data:
        ls.append(deStr(i))
    data = sorted(ls)
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
                # Dict[str(count) + " " + header] = seq
                Dict[count]["seq"] = seq
                Dict[count]["header"] = header
                count += 1
                header = i[1:]
                seq = ''
            else:
                header = i[1:]
                seq = ''
        else:
            seq += i
    Dict[count]["seq"] = seq
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
    Dict[header]["seq"] = seq
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


def main():
    if args.db_type == "nr" or args.db_type == "refseq":
        os.system("mkdir " + args.output + "/blast")
        print("made file: " + args.output + "/blast")
        print("blasting against the NCBI database...")
        os.system("blastp -query " + args.genes + " -db " + args.db + " -num_threads " + str(args.t) + " -out " + cwd +
                  "/" + args.output + "/blast/blastResults -evalue " + str(args.e) + " -outfmt 6")
        print("blast finished\n\nreading in the results\n")
        blast = open(cwd + "/" + args.output + "/blast/blastResults", "r")
        AccessionDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        bDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in blast:
            bList = i.rstrip().split("\t")
            accession = bList[1]
            bDict[accession] = bList
            AccessionDict[deStr(accession)] = accession

        print("identifying gene cassettes...")
        x = cluster(bDict.keys(), maxgap=args.g)
        outfile = open(cwd + "/" + args.output + "/syntenies.csv", "w")
        outfile.write(
            "gene" + "," + "qseqid" + "," + "sseqid" + "," + "pident" + "," + "length" + "," + "mismatch" + ","
            + "gapopen" + "," + "qstart" + "," + "qend" + "," + "sstart" + "," + "send" + "," + "evalue" +
            "," + "bitscore" + "\n")
        for i in x:
            if len(i) >= int(args.n):
                print("found cassette: " + str(i))
                print("identifying genes from accession numbers")
                for accession in i:
                    try:
                        fp = urllib.request.urlopen(
                            "https://www.ncbi.nlm.nih.gov/protein/%s" % AccessionDict[accession], context=gcontext)
                        mybytes = fp.read()
                        mystr = mybytes.decode("utf8")
                        fp.close()
                        lines = mystr.split('\n')
                        for j in lines:
                            if re.findall(r'<h1>', j):
                                # try:
                                    j1 = (j.split(">")[1])
                                    j2 = (j1.split("<")[0])
                                    j2 = remove(j2, [","])
                                    # org = (re.search(r'((\[)(.*)(\]))', j2).group(1))
                                    # prot = (j2.split("[")[0])
                                    # outfile.write(prot + "," + org + ",")
                                    outfile.write(j2 + ",")
                                    accession = AccessionDict[accession]
                                    accession = str(accession)
                                    for k in bDict[accession]:
                                        outfile.write(k + ",")
                                    outfile.write("\n")
                                # except IndexError:
                                #     outfile.write(j2 + "," + "" + ",")
                                #     for col in bDict[accession]:
                                #         outfile.write(col + ",")
                                #     outfile.write("\n")
                            outfile.write("\n")
                    except HTTPError:
                        outfile.write("Not found in NCBI" + ",")
                        accession = AccessionDict[accession]
                        accession = str(accession)
                        for k in bDict[accession]:
                            outfile.write(k + ",")
                        outfile.write("\n")
    # ================================================================================================================
    # If there is a local directory provided
    else:
        outfile = open(cwd + "/" + args.output + "/syntenies.csv", "w")
        outfile.write(
            "gene" + "," + "qseqid" + "," + "sseqid_assigned" + "," + "pident" + "," + "length" + "," + "mismatch" + ","
            + "gapopen" + "," + "qstart" + "," + "qend" + "," + "sstart" + "," + "send" + "," + "evalue" +
            "," + "bitscore" + "\n")
        os.system("mkdir " + args.output + "/blast")
        FastaDir = os.listdir(args.db)
        for genome in FastaDir:
            print(genome)
            Genome = open(args.db + "/" + genome, "r")
            Genome = fasta(Genome)
            out = open(args.db + "/" + genome + ".numbered", "w")
            for key in (Genome.keys()):
                out.write(">" + str(key) + "\n")
                out.write(Genome[key]['seq'] + "\n")
            out.close()
            os.system(
                "makeblastdb -dbtype prot -in " + args.db + "/" + genome + ".numbered -out " + args.db + "/" + genome + ".numbered -logfile blastLog")
            os.system(
                "blastp -query " + args.genes + " -db " + args.db + "/" + genome + ".numbered -num_threads " + str(
                    args.t) + " -out " +
                cwd + "/" + args.output + "/blast/" + genome + ".blastResults -evalue " + str(args.e) + " -outfmt 6")
            blast = open(cwd + "/" + args.output + "/blast/" + genome + ".blastResults", "r")
            bDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in blast:
                bList = i.rstrip().split("\t")
                num = bList[1]
                bDict[num] = bList
            x = cluster(bDict.keys(), maxgap=args.g)
            for i in x:
                if len(i) >= int(args.n):
                    for num in i:
                        num = int(num)
                        outfile.write(remove(Genome[num]["header"], [","]) + ",")
                        num = str(num)
                        for k in bDict[num]:
                            outfile.write(k + ",")
                        outfile.write("\n")
                    outfile.write("\n")
            outfile.write("\n")

            os.system("rm " + args.db + "/" + genome + ".numbered.phr")
            os.system("rm " + args.db + "/" + genome + ".numbered.pin")
            os.system("rm " + args.db + "/" + genome + ".numbered.psq")
            os.system("rm " + args.db + "/" + genome + ".numbered")
            os.system("rm blastLog")
        outfile.close()
        os.system("rm blastLog.perf")


cwd = os.getcwd()
os.system("mkdir " + cwd + "/" + args.output)
main()
