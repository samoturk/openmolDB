import pybel
import openbabel
import csv

from openmoldbmolecules.models import Molecule

def make_db_full(data):
    namec = ""
    name2c = ""
    suppc = ""
    suppidc = ""
    storc = ""
    amountc = ""
    unitc = ""
    casc = ""
    supp2c = ""
    suppID2c = ""
    commc = ""
    classc = ""
    counter = 0
    head = data.next()
    for c in head:
        #Get the indexes of columns
        cl = c.lower()
        if cl == "name":
            namec = counter
        elif cl == "name2":
            name2c = counter
        elif cl == "supplier":
            suppc = counter
        elif cl == "supplierid":
            suppidc = counter
        elif cl == "storageid":
            storc = counter
        elif cl == "amount":
            amountc = counter
        elif cl == "unit":
            unitc = counter
        elif cl == "cas":
            casc = counter
        elif cl == "smiles":
            smilesc = counter
        elif cl == "supplier2":
            supp2c = counter
        elif cl == "supplierid2":
            suppID2c = counter
        elif cl == "comment":
            commc = counter
        elif cl == "class":
            classc = counter
        counter += 1

    if type(smilesc) == int or type(namec) == int:
        # has to have either smiles or name column
        print "Following columns were found and will be imported:"
        if type(namec) == int:
            print head[namec] + " ",
            
        if type(name2c) == int:
            print head[name2c] + " ",
            
        if type(suppc) == int:
            print head[suppc] + " ",
            
        if type(suppidc) == int:
            print head[suppidc] + " ",
            
        if type(storc) == int:
            print head[storc] + " ",
            
        if type(amountc) == int:
            print head[amountc] + " ",
            
        if type(unitc) == int:
            print head[unitc] + " ",
            
        if type(casc) == int:
            print head[casc] + " ",
            
        if type(smilesc) == int:
            print head[smilesc] + " ",
        
        if type(supp2c) == int:
            print head[supp2c] + " ",
        
        if type(suppID2c) == int:
            print head[suppID2c] + " ",
        
        if type(commc) == int:
            print head[commc] + " ",
        
        if type(classc) == int:
            print head[classc] + " ",
        
        print ""
            
        userinput =  raw_input("Is this ok? yes/no: ")
        #Check with user if everything looks OK
        
        if userinput == "yes":
            for line in data:
                #do some datachecks and encode in ACSII since some databases have problems with UTF-8
                if type(namec) == int:
                    name = line[namec].decode("windows-1252").encode('utf-8','ignore')
                else:
                    name = ""
                
                if type(name2c) == int:
                    name2 = line[name2c].decode("windows-1252").encode('utf-8','ignore')
                else:
                    name2 = ""
            
                if type(suppc) == int:
                    supp = line[suppc].decode("windows-1252").encode('utf-8','ignore')
                else:
                    supp = ""
            
                if type(suppidc) == int:
                    suppid = line[suppidc].decode("windows-1252").encode('utf-8','ignore')
                else:
                    suppid = ""
            
                if type(storc) == int:
                    storageid = line[storc].decode("windows-1252").encode('utf-8','ignore')
                else:
                    storageid = ""
            
                if type(amountc) == int:
                    amount = line[amountc]
                else:
                    amount = ""
            
                if type(unitc) == int:
                    unit = line[unitc].decode("windows-1252").encode('utf-8','ignore')
                else:
                    unit = ""
            
                if type(casc) == int:
                    cas = line[casc].decode("windows-1252").encode('utf-8','ignore')
                else:
                    cas = ""
            
                if type(smilesc) == int:
                    smiles = line[smilesc].decode("windows-1252").encode('utf-8','ignore')
                else:
                    smiles = ""
        
                if type(supp2c) == int:
                    supp2 = line[supp2c].decode("windows-1252").encode('utf-8','ignore')
                else:
                    supp2 = ""
        
                if type(suppID2c) == int:
                    suppid2 = line[suppID2c].decode("windows-1252").encode('utf-8','ignore')
                else:
                    suppid2 = ""
        
                if type(commc) == int:
                    comm = line[commc].decode("windows-1252").encode('utf-8','ignore')
                else:
                    comm = ""
                
                if type(classc) == int:
                    typeofc = line[classc].decode("windows-1252").encode('utf-8','ignore')
                    print typeofc
                else:
                    typeofc = "chemical"
            
                try:
                    amount = float(line[amountc])
                    unit = line[unitc]
                except:
                    amount = 0
                    unit = "X"

                try:
                    mol = pybel.readstring("smi",smiles)
                    descs = mol.calcdesc()
                    #generate 2D coordinates, needs openbabel
                    obConversion = openbabel.OBConversion()
                    obConversion.SetInAndOutFormats("smi", "mdl")
                    obmol = openbabel.OBMol()
                    obConversion.ReadString(obmol, smiles)
                    gen2d = openbabel.OBOp.FindType("gen2d")
                    gen2d.Do(obmol)
                    MDL = obConversion.WriteString(obmol)
                    outMDL = MDL.replace("\n", r"\n")
                    CMW=descs["MW"]
                    CHN=mol.formula
                    HBA=descs["HBA1"]
                    HBD=descs["HBD"]
                    logP=descs["logP"]
                    tpsa=descs["TPSA"]
                    #Get number  of rotatable bonds
                    smarts = pybel.Smarts("[!$([NH]!@C(=O))&!D1&!$(*#*)]\&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]")
                    rb = smarts.findall(mol)
                    nrb = len(rb)
                    #Get fingerprint and molecular complexity
                    fprint = mol.calcfp()
                    bitson = fprint.bits
                    nbitson = len(bitson)
                    print name
                    m = Molecule(name=name,SMILES=smiles, name2=name2, supp1=supp, suppID1=suppid, CMW=descs["MW"], CHN=CHN, HBA=HBA, HBD=HBD, logP=logP, tpsa=tpsa, amount=amount, unit=unit, CAS=cas, storageID=storageid, molfile=outMDL, nrb=nrb, fingerprint=bitson, complexity=nbitson, supp2=supp2, suppID2=suppid2 , comment=comm, typeOfCompound=typeofc)
                    m.save()
                except:
                    m = Molecule(name=name,SMILES=smiles, name2=name2, supp1=supp, suppID1=suppid, amount=amount, unit=unit, CAS=cas, storageID=storageid, supp2=supp2, suppID2=suppid2 , comment=comm, typeOfCompound=typeofc)
                    m.save()
            #Save data to database
        else:
            print "Exiting, no changes were made..."
            return False
    else:
        print "No valid columns were found in the table"

if __name__ == '__main__':	
    import sys
    datacard = sys.argv[1]
    data = csv.reader(open(datacard), delimiter='\t')
    make_db_full(data)
