import pybel
import openbabel
import csv

from openmoldbmolecules.models import Molecule
from openmoldbmolecules.filter_pains import detect_pains

def addsingle(name, altname, supplier, supplierID, storage, storageID, unit, amount, cas, smiles, comment, molclass, platebarcode, samplebarcode, randomstring):
    #do some datachecks and encode in ACSII since some databases have problems with UTF-8
    name = name.decode("windows-1252").encode('utf-8','ignore')
    altname = altname.decode("windows-1252").encode('utf-8','ignore')
    supplier = supplier.decode("windows-1252").encode('utf-8','ignore')
    supplierID = supplierID.decode("windows-1252").encode('utf-8','ignore')
    storage = storage.decode("windows-1252").encode('utf-8','ignore')
    storageID = storageID.decode("windows-1252").encode('utf-8','ignore')
    unit = unit.decode("windows-1252").encode('utf-8','ignore')
    cas = cas.decode("windows-1252").encode('utf-8','ignore')
    smiles = smiles.decode("windows-1252").encode('utf-8','ignore')
    comment = comment.decode("windows-1252").encode('utf-8','ignore')
    molclass = molclass.decode("windows-1252").encode('utf-8','ignore')
    platebarcode = platebarcode.decode("windows-1252").encode('utf-8','ignore')
    samplebarcode = samplebarcode.decode("windows-1252").encode('utf-8','ignore')
    # Make sure amount is a number
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
        smarts = pybel.Smarts(r"[!$([NH]!@C(=O))&!D1&!$(*#*)]\&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]")
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
        else:
            fsp3 = ""
        #Get fingerprint and molecular complexity
        fprint = mol.calcfp()
        bitson = fprint.bits
        nbitson = len(bitson)
        if 'hts' in molclass.lower() or 'compound' in molclass.lower():
            #print "hts"
            pains = detect_pains(mol)
        else:
            pains = 'Not checked'
        m = Molecule(name=name,SMILES=smiles, altname=altname, supplier=supplier, supplierID=supplierID, 
        CMW=descs["MW"], CHN=CHN, HBA=HBA, HBD=HBD, logP=logP, tpsa=tpsa, amount=amount, unit=unit, 
        CAS=cas, storage=storage, storageID=storageID, molfile=outMDL, nrb=nrb, fingerprint=bitson, complexity=nbitson, 
        comment=comment, molclass=molclass, fsp3=fsp3, pains=pains, platebarcode=platebarcode, 
        samplebarcode=samplebarcode, randomstring=randomstring)
        m.save()
    except:
        # OpenBabel failed, no properties, etc..
        m = Molecule(name=name,SMILES=smiles, altname=altname, supplier=supplier, supplierID=supplierID, 
        amount=amount, unit=unit, CAS=cas, storage=storage, storageID=storageID, comment=comment, molclass=molclass, 
        platebarcode=platebarcode, samplebarcode=samplebarcode, randomstring=randomstring)
        m.save()
        #Save data to database

def addtable(csvdata, cmdline, userinput, randomstring):
    """
    accepts csv table (from csv module)
    checks first row for colums (name, altname, smiles, etc)
    adds each row to the database
    """
    namec = ""
    altnamec = ""
    suppc = ""
    suppidc = ""
    storc = ""
    storidc = ""
    amountc = ""
    unitc = ""
    casc = ""
    altsupplierc = ""
    altsupplierIDc = ""
    commc = ""
    classc = ""
    platec = ""
    samplec = ""
    counter = 0
    head = csvdata.next()
    print "Found columns:"
    for c in head:
        #Get the indexes of columns
        cl = c.lower()
        if cl == "name":
            namec = counter
            print "name ",
        elif cl == "altname":
            altnamec = counter
            print "altname ",
        elif cl == "supplier":
            suppc = counter
            print "supplier ",
        elif cl == "supplierid":
            suppidc = counter
            print "supplierid ",
        elif cl == "storage":
            storc = counter
            print "storage ",
        elif cl == "storageid":
            storidc = counter
            print "storageid ",
        elif cl == "amount":
            amountc = counter
            print "amount ",
        elif cl == "unit":
            unitc = counter
            print "unit ",
        elif cl == "cas":
            casc = counter
            print "cas ",
        elif cl == "smiles":
            smilesc = counter
            print "smiles ",
        elif cl == "supplier2":
            altsupplierc = counter
            print "supplier2 ",
        elif cl == "supplierid2":
            altsupplierIDc = counter
            print "supplierid2 ",
        elif cl == "comment":
            commc = counter
            print "comment ",
        elif cl == "molclass":
            classc = counter
            print "molclass ",
        elif cl == "platebarcode":
            platec = counter
            print "platebarcode ",
        elif cl == "samplebarcode":
            samplec = counter
            print "samplebarcode ",
        counter += 1

    if type(smilesc) == int or type(namec) == int:
        # has to have either smiles or name column
        if cmdline:
            print ""
            userinput =  raw_input("Is this ok? yes/no: ")
        #Check with user if everything looks OK
        
        if userinput == "yes":
            for line in csvdata:
                if type(namec) == int:
                    name = line[namec]
                else:
                    name = ""
                
                if type(altnamec) == int:
                    altname = line[altnamec]
                else:
                    altname = ""
            
                if type(suppc) == int:
                    supplier = line[suppc]
                else:
                    supplier = ""
            
                if type(suppidc) == int:
                    supplierID = line[suppidc]
                else:
                    supplierID = ""
            
                if type(storc) == int:
                    storage = line[storc]
                else:
                    storage = ""
                
                if type(storidc) == int:
                    storageID = line[storidc]
                else:
                    storageID = ""
            
                if type(amountc) == int:
                    amount = line[amountc]
                else:
                    amount = ""
            
                if type(unitc) == int:
                    unit = line[unitc]
                else:
                    unit = ""
            
                if type(casc) == int:
                    cas = line[casc]
                else:
                    cas = ""
            
                if type(smilesc) == int:
                    smiles = line[smilesc]
                else:
                    smiles = ""
        
                if type(altsupplierc) == int:
                    altsupplier = line[altsupplierc]
                else:
                    altsupplier = ""
        
                if type(altsupplierIDc) == int:
                    altsupplierID = line[altsupplierIDc]
                else:
                    altsupplierID = ""
        
                if type(commc) == int:
                    comment = line[commc]
                else:
                    comment = ""
                
                if type(classc) == int:
                    molclass = line[classc]
                else:
                    molclass = ""
                
                if type(platec) == int:
                    platebarcode = line[platec]
                else:
                    platebarcode = ""
                
                if type(samplec) == int:
                    samplebarcode = line[samplec]
                else:
                    samplebarcode = ""
            
                try:
                    amount = float(line[amountc])
                    unit = line[unitc]
                except:
                    amount = 0
                    unit = "X"
                
                addsingle(name, altname, supplier, supplierID, storage, storageID, unit, amount, cas, smiles, comment, molclass, platebarcode, samplebarcode, randomstring)
            
            return "OK"
        else:
            return "No mols were added. Exiting!"
    else:
        return "No valid columns were found"