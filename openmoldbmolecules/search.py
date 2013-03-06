from openmoldbmolecules.models import Molecule
import pybel
from sets import Set

def smiles_search(mollist, smiles, tanimoto):
    ret = []
    tanimoto = float(tanimoto)
    query = pybel.readstring("smi", smiles)
    fpquery = query.calcfp()
    for mol in mollist.all():
        try:
            smiles = pybel.readstring("smi", str(mol.SMILES))
            fpsmiles = smiles.calcfp()
            tan = fpquery | fpsmiles
            #print tan
            if tan > tanimoto:
                ret.append(mol)
        except:
            pass
    return ret

def fast_fp_search(mollist, smiles, tanimoto):
    #Should be faster, since it uses precalculated fingerprints
    ret = []
    tanimoto = float(tanimoto)
    query = pybel.readstring("smi", smiles)
    fpquery = query.calcfp()
    fpquery = fpquery.bits
    for mol in mollist.all():
        try:
            fp = eval(str(mol.fingerprint))
            a = Set(fpquery)
            b = Set(fp)
            un = float(len(a.union(b)))
            inx = float(len(a.intersection(b)))
            tan = inx/un
            #print tan
            if tan > tanimoto:
                ret.append(mol)
        except:
            pass
    return ret

def smarts_search(mollist, smarts):
    ret = []
    query = pybel.Smarts(smarts)
    print query
    for mol in mollist.all():
        try:
            smiles = pybel.readstring("smi", str(mol.SMILES))
            if query.findall(smiles):
                ret.append(mol)
        except:
            pass
    return ret

def properties_search(mollist, minmw, maxmw, minlogp, maxlogp, minhba, maxhba, minhbd, maxhbd):
    ret = []
    for mol in mollist:
        if all([mol.CMW >= minmw, mol.CMW <= maxmw, minlogp <= mol.logP,  mol.logP <= maxlogp, minhba <= mol.HBA,  mol.HBA <= maxhba , minhbd <= mol.HBD,  mol.HBD <= maxhbd]):
            ret.append(mol)
    return ret

def id_search(cas, name, storageid, supplierid, supplier, chn, molclass):
    if cas:
        #CAS should be unique, so no need for further filtering
        mols = Molecule.objects.filter(CAS__contains=cas)
        return mols
    else:
        mols = Molecule.objects.all()
        if name:
            mols = mols.filter(name__icontains=name)
        if storageid:
            mols = mols.filter(storageID__icontains=storageid)
        if supplierid:
            mols = mols.filter(supplierID__icontains=supplierid)
        if supplier:
            mols = mols.filter(supplier__icontains=supplier)
        if chn:
            mols = mols.filter(CHN__icontains=chn)
        if molclass:
            mols = mols.filter(molclass__icontains=molclass)
        return mols
    