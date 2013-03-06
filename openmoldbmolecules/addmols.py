import pybel
import openbabel
import csv

from openmoldbmolecules.models import Molecule

def addsingle(name, altname, supplier, supplierID, storageID, unit, amount, cas, smiles, comment, molclass, randomstring):
    #do some datachecks and encode in ACSII since some databases have problems with UTF-8
    name = name.decode("windows-1252").encode('utf-8','ignore')
    altname = altname.decode("windows-1252").encode('utf-8','ignore')
    supplier = supplier.decode("windows-1252").encode('utf-8','ignore')
    supplierID = supplierID.decode("windows-1252").encode('utf-8','ignore')
    storageID = storageID.decode("windows-1252").encode('utf-8','ignore')
    unit = unit.decode("windows-1252").encode('utf-8','ignore')
    cas = cas.decode("windows-1252").encode('utf-8','ignore')
    smiles = smiles.decode("windows-1252").encode('utf-8','ignore')
    comment = comment.decode("windows-1252").encode('utf-8','ignore')
    molclass = molclass.decode("windows-1252").encode('utf-8','ignore')
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
        else:
            fsp3 = ""
        #Get fingerprint and molecular complexity
        fprint = mol.calcfp()
        bitson = fprint.bits
        nbitson = len(bitson)
        m = Molecule(name=name,SMILES=smiles, altname=altname, supplier=supplier, supplierID=supplierID, CMW=descs["MW"], CHN=CHN, HBA=HBA, HBD=HBD, logP=logP, tpsa=tpsa, amount=amount, unit=unit, CAS=cas, storageID=storageID, molfile=outMDL, nrb=nrb, fingerprint=bitson, complexity=nbitson, comment=comment, molclass=molclass, fsp3=fsp3, randomstring=randomstring)
        m.save()
    except:
        # OpenBabel failed, no properties, etc..
        m = Molecule(name=name,SMILES=smiles, altname=altname, supplier=supplier, supplierID=supplierID, amount=amount, unit=unit, CAS=cas, storageID=storageID, comment=comment, molclass=molclass, randomstring=randomstring)
        m.save()
        #Save data to database

