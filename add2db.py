import csv
from openmoldbmolecules.addmols import addtable

def add2db(csvfile):
    cmdline = True
    userinput="yes"
    randomstring = ""
    csvdata = csv.reader(open(csvfile), dialect="excel")
    print addtable(csvdata, cmdline, userinput, randomstring)