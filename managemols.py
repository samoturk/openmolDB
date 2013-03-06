import pybel
import openbabel
import csv

from openmoldbmolecules.models import Molecule

def add_mols(data):
    namec = ""
    altnamec = ""
    suppc = ""
    suppidc = ""
    storc = ""
    amountc = ""
    unitc = ""
    casc = ""
    altsupplierc = ""
    altsupplierIDc = ""
    commc = ""
    classc = ""
    counter = 0
    head = data.next()
    for c in head:
        #Get the indexes of columns
        cl = c.lower()
        if cl == "name":
            namec = counter
        elif cl == "altname":
            altnamec = counter
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
            altsupplierc = counter
        elif cl == "supplierid2":
            altsupplierIDc = counter
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
            
        if type(altnamec) == int:
            print head[altnamec] + " ",
            
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
        
        if type(altsupplierc) == int:
            print head[altsupplierc] + " ",
        
        if type(altsupplierIDc) == int:
            print head[altsupplierIDc] + " ",
        
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
                
                if type(altnamec) == int:
                    altname = line[altnamec].decode("windows-1252").encode('utf-8','ignore')
                else:
                    altname = ""
            
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
        
                if type(altsupplierc) == int:
                    altsupplier = line[altsupplierc].decode("windows-1252").encode('utf-8','ignore')
                else:
                    altsupplier = ""
        
                if type(altsupplierIDc) == int:
                    altsupplierID = line[altsupplierIDc].decode("windows-1252").encode('utf-8','ignore')
                else:
                    altsupplierID = ""
        
                if type(commc) == int:
                    comm = line[commc].decode("windows-1252").encode('utf-8','ignore')
                else:
                    comm = ""
                
                if type(classc) == int:
                    molclass = line[classc].decode("windows-1252").encode('utf-8','ignore')
                else:
                    molclass = ""
            
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
                    #Calculate Fsp3
                    sp3c = pybel.Smarts("[CX4]")
                    nsp3c = sp3c.findall(mol)
                    nsp3c = float(len(nsp3c))
                    allc =  pybel.Smarts("[#6]")
                    nallc = allc.findall(mol)
                    nallc = float(len(nallc))
                    if nallc > 0:
                        fsp3 = nsp3c/nallc
                        print fsp3
                    else:
                        fsp3 = ""
                    #Get fingerprint and molecular complexity
                    fprint = mol.calcfp()
                    bitson = fprint.bits
                    nbitson = len(bitson)
                    print name
                    m = Molecule(name=name,SMILES=smiles, altname=altname, supplier=supp, supplierID=suppid, CMW=descs["MW"], CHN=CHN, HBA=HBA, HBD=HBD, logP=logP, tpsa=tpsa, amount=amount, unit=unit, CAS=cas, storageID=storageid, molfile=outMDL, nrb=nrb, fingerprint=bitson, complexity=nbitson, altsupplier=altsupplier, altsupplierID=altsupplierID , comment=comm, molclass=molclass, fsp3=fsp3)
                    m.save()
                except:
                    m = Molecule(name=name,SMILES=smiles, altname=altname, supplier=supp, supplierID=suppid, amount=amount, unit=unit, CAS=cas, storageID=storageid, altsupplier=altsupplier, altsupplierID=altsupplierID , comment=comm, molclass=molclass)
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
    data = csv.reader(open(datacard), dialect='excel')
    add_mols(data)
