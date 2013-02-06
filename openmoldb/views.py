#!/usr/bin/python/
from django.shortcuts import render, render_to_response, get_object_or_404
from openmoldbmolecules.models import Molecule
from openmoldbmolecules.search import smiles_search, properties_search, id_search, smarts_search, fast_fp_search
from openmoldbmolecules.forms import SearchForm
from openmoldbmolecules.pagination import get_page_wo_page
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from json import dumps
from numpy import histogram, array
import pybel
#from django.http import HttpResponse 

servername = "http://127.0.0.1:8000"

def home(request):
    #home page
    return render(request, 'index.html', {'servername':servername})

def molecule(request, idnum):
    #Views for individual molecules
    try:
        offset = int(idnum)
    except ValueError:
        raise Http404()
    mol = get_object_or_404(Molecule, pk=idnum)
    quantity = str(mol.amount) + " " + str(mol.unit)
    return render(request, 'molecule.html', {'servername':servername, 'mol':mol, 'molname':mol.name, 'smiles':mol.SMILES, 'cas':mol.CAS, 'name2':mol.name2, 'molfile':mol.molfile, 'mw':mol.CMW, 'chn':mol.CHN, 'hba':mol.HBA, 'hbd':mol.HBD, 'logp':mol.logP, 'tpsa':mol.tpsa, 'supp':mol.supp1, 'suppid':mol.suppID1, 'storage':mol.storageID, 'comm':mol.comment, 'quantity':quantity, 'added':mol.added})

def browse(request):
    mollist = Molecule.objects.all()
    paginator = Paginator(mollist, 15)
    page = request.GET.get('page') #get value of query
    nofmols = len(mollist)
    try:
        mollist = paginator.page(page)
    except PageNotAnInteger:
        # If page is not an integer, deliver first page.
        mollist = paginator.page(1)
    except EmptyPage:
        # If page is out of range (e.g. 9999), deliver last page of results.
        mollist = paginator.page(paginator.num_pages)
    return render(request, 'browse.html', {'servername':servername, 'mollist':mollist, 'currpage':page, 'ofmols':nofmols, 'currentpath':'?'})

def search(request):
    error = []
    form = SearchForm(request.GET)
    form.is_valid()
    # Basic form validation as defined in forms.py
    if form.is_valid() == False:
        for item in form.errors:
            error.append(form[item].errors)
        return render(request, 'search.html', {'servername':servername, 'error': error, 'form':form})
    moltext = form.cleaned_data['moltext']
    minmw = form.cleaned_data['minmw']
    maxmw = form.cleaned_data['maxmw']
    minlogp = form.cleaned_data['minlogp']
    maxlogp = form.cleaned_data['maxlogp']
    minhba = form.cleaned_data['minhba']
    maxhba = form.cleaned_data['maxhba']
    minhbd = form.cleaned_data['minhbd']
    maxhbd = form.cleaned_data['maxhbd']
    chn = form.cleaned_data['chn']
    name = form.cleaned_data['name']
    storageid = form.cleaned_data['storageid']
    supplierid = form.cleaned_data['supplierid']
    supplier = form.cleaned_data['supplier']
    cas = form.cleaned_data['cas']
    if 'similarity' in request.GET:
        #similartiy search
        try:
            #check if smiles or mol from editor is valid and run search
            if moltext:
                try:
                    moltext = pybel.readstring("mol", str(moltext))
                    smiles = moltext.write("smi")
                except:
                    smiles = str(form.cleaned_data['smiles'])
            query = pybel.readstring("smi", smiles)
            tanimoto = form.cleaned_data['tanimoto']
            mollist = id_search(cas, name, storageid, supplierid, supplier, chn)
            mollistsmiles = fast_fp_search(mollist, smiles, tanimoto)
            mollist = properties_search(mollistsmiles, minmw, maxmw, minlogp, maxlogp, minhba, maxhba, minhbd, maxhbd)
            paginator = Paginator(mollist, 15)
            page = request.GET.get('page') #get value of query
            nofmols = len(mollist)
            currentpath = request.get_full_path()
            currentpath = get_page_wo_page(currentpath)
            try:
                mollist = paginator.page(page)
            except PageNotAnInteger:
                # If page is not an integer, deliver first page.
                mollist = paginator.page(1)
            except EmptyPage:
                # If page is out of range (e.g. 9999), deliver last page of results.
                mollist = paginator.page(paginator.num_pages)
            return render(request, 'browse.html', {'servername':servername, 'mollist':mollist,'currpage':page, 'ofmols':nofmols, 'currentpath':currentpath})
        except:
            #if smiles is not valid, return error
            error.append('SMILES code or tanimoto index is not valid!')
            return render(request, 'search.html', {'servername':servername, 'error': error, 'form':form})
    elif 'substructure' in request.GET:
        #substructure search
        try:
            #check if smiles/smarts or mol from editor is valid
            if moltext:
                try:
                    moltext = pybel.readstring("mol", str(moltext))
                    smarts = moltext.write("smi")
                except:
                    smarts = str(form.cleaned_data['smiles'])
            query = pybel.Smarts(smarts)
            mollist = id_search(cas, name, storageid, supplierid, supplier, chn)
            mollistsmart = smarts_search(mollist, smarts)
            mollist = properties_search(mollistsmart, minmw, maxmw, minlogp, maxlogp, minhba, maxhba, minhbd, maxhbd)
            paginator = Paginator(mollist, 15)
            page = request.GET.get('page') #get value of query
            nofmols = len(mollist)
            currentpath = request.get_full_path()
            currentpath = get_page_wo_page(currentpath)
            try:
                mollist = paginator.page(page)
            except PageNotAnInteger:
                # If page is not an integer, deliver first page.
                mollist = paginator.page(1)
            except EmptyPage:
                # If page is out of range (e.g. 9999), deliver last page of results.
                mollist = paginator.page(paginator.num_pages)
            return render(request, 'browse.html', {'servername':servername, 'mollist':mollist,'currpage':page, 'ofmols':nofmols, 'currentpath':currentpath})
        except:
            #if smiles is not valid, return error
            error.append('SMILES/SMARTS code is not valid!')
            return render(request, 'search.html', {'servername':servername, 'error': error, 'form':form})
    elif 'properties' in request.GET:
        #properties and IDs search
        mollist = id_search(cas, name, storageid, supplierid, supplier, chn)
        if cas == '':
            mollist = properties_search(mollist, minmw, maxmw, minlogp, maxlogp, minhba, maxhba, minhbd, maxhbd)
        
        paginator = Paginator(mollist, 15)
        page = request.GET.get('page') #get value of query
        nofmols = len(mollist)
        currentpath = request.get_full_path()
        currentpath = get_page_wo_page(currentpath)
        
        try:
            mollist = paginator.page(page)
        except PageNotAnInteger:
            # If page is not an integer, deliver first page.
            mollist = paginator.page(1)
        except EmptyPage:
            # If page is out of range (e.g. 9999), deliver last page of results.
            mollist = paginator.page(paginator.num_pages)
        return render(request, 'browse.html', {'servername':servername, 'mollist':mollist,'currpage':page, 'ofmols':nofmols, 'currentpath':currentpath})
    
    else:
        form = SearchForm()
    return render(request, 'search.html', {'servername':servername, 'error': error, 'form':form})
    
def statistics(request):
    logp = []
    mw = []
    hba = []
    hbd = []
    tpsa = []
    names = []
    fsp3 = []
    nrb = []
    complexity = []
    #c = 0 # uncomment for testing, it stops after 100
    for mol in Molecule.objects.all():
        try:
            float(mol.logP)
            float(mol.CMW)
            float(mol.HBA)
            float(mol.HBD)
            float(mol.tpsa)
            float(mol.fsp3)
            float(mol.nrb)
            float(mol.complexity)
            if mol.CMW <= 800:
                logp.append(mol.logP)
                mw.append(mol.CMW)
                hba.append(mol.HBA)
                hbd.append(mol.HBD)
                tpsa.append(mol.tpsa)
                fsp3.append(mol.fsp3*100)
                nrb.append(mol.nrb)
                complexity.append(mol.complexity)
                if mol.name == "":
                    names.append("X")
                else:
                    names.append(mol.name)
            #logpmw +=  "[" + str(mol.logP) + ", " + str(mol.CMW)+ ", '" + str(molname) + "'], "
            #c += 1 # uncomment for testing
        except:
            pass
        #if c == 100: # uncomment for testing
        #    break # uncomment for testing
    #logpmw = logpmw[:-2]
    logpmw = dumps(zip(logp, mw))
    hbamw = dumps(zip(hba, mw))
    hbdmw = dumps(zip(hbd, mw))
    tpsamw = dumps(zip(tpsa, mw))
    #histograms
    mwhist, mwedge = histogram(mw, bins=20)
    mwhist = dumps(zip(array(mwedge[:-1]).tolist(), array(mwhist).tolist()))
    fsp3hist, fsp3edge = histogram(fsp3, bins=20)
    fsp3hist = dumps(zip(array(fsp3edge[:-1]).tolist(), array(fsp3hist).tolist()))
    nrbhist, nrbedge = histogram(nrb, bins=35)
    nrbhist = dumps(zip(array(nrbedge[:-1]).tolist(), array(nrbhist).tolist()))
    complexityhist, complexityedge = histogram(complexity, bins=20)
    complexityhist = dumps(zip(array(complexityedge[:-1]).tolist(), array(complexityhist).tolist()))
    logphist, logpedge = histogram(logp, bins=20)
    logphist = dumps(zip(array(logpedge[:-1]).tolist(), array(logphist).tolist()))
    hbahist, hbaedge = histogram(hba, bins=20)
    hbahist = dumps(zip(array(hbaedge[:-1]).tolist(), array(hbahist).tolist()))
    hbdhist, hbdedge = histogram(hbd, bins=15)
    hbdhist = dumps(zip(array(hbdedge[:-1]).tolist(), array(hbdhist).tolist()))
    tpsahist, tpsaedge = histogram(tpsa, bins=15)
    tpsahist = dumps(zip(array(tpsaedge[:-1]).tolist(), array(tpsahist).tolist()))
    return render(request, 'statistics.html', {'servername':servername, 'logpmw':logpmw, 'hbamw':hbamw, 'hbdmw':hbdmw, 'tpsamw':tpsamw, 'mwhist':mwhist, 'fsp3hist':fsp3hist, 'nrbhist':nrbhist, 'complexityhist':complexityhist, 'logphist':logphist, 'hbahist':hbahist, 'hbdhist':hbdhist, 'tpsahist':tpsahist})

