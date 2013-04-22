#!/usr/bin/python/
from django.shortcuts import render, render_to_response, get_object_or_404, HttpResponseRedirect
from openmoldbmolecules.models import Molecule
from openmoldbmolecules.search import smiles_search, properties_search, id_search, smarts_search, fast_fp_search
from openmoldbmolecules.forms import SearchForm, SubmitSingle, UploadFileForm
from openmoldbmolecules.pagination import get_page_wo_page
from openmoldbmolecules.addmols import addsingle, addtable
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.contrib import auth
from json import dumps
from numpy import histogram, array
import pybel
import string
#import random
import csv
#from django.http import HttpResponse 

#Some settings
servername = "http://127.0.0.1:8000" #Change to your server name

logintoupload = True # Require login to upload
uploadsingle = True # Option to upload single molecule
uploadtable = True # Option to upload csv table with mols
checkpains = True # Option to check pains if mol is compound or hts

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
    return render(request, 'molecule.html', {'servername':servername, 'mol':mol, 'molname':mol.name, 'smiles':mol.SMILES, 'cas':mol.CAS, 'altname':mol.altname, 'molfile':mol.molfile, 'mw':mol.CMW, 'chn':mol.CHN, 'hba':mol.HBA, 'hbd':mol.HBD, 'logp':mol.logP, 'tpsa':mol.tpsa, 'supp':mol.supplier, 'suppid':mol.supplierID, 'storage':mol.storageID, 'comm':mol.comment, 'quantity':quantity, 'added':mol.added})

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
    mintpsa = form.cleaned_data['mintpsa']
    maxtpsa = form.cleaned_data['maxtpsa']
    minfsp3 = form.cleaned_data['minfsp3']
    maxfsp3 = form.cleaned_data['maxfsp3']
    chn = form.cleaned_data['chn']
    name = form.cleaned_data['name']
    storageid = form.cleaned_data['storageid']
    supplierid = form.cleaned_data['supplierid']
    supplier = form.cleaned_data['supplier']
    cas = form.cleaned_data['cas']
    molclass = form.cleaned_data['molclass']
    platebarcode = form.cleaned_data['platebarcode']
    samplebarcode = form.cleaned_data['samplebarcode']
    storage = form.cleaned_data['storage']
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
            mollist = id_search(cas, name, storageid, supplierid, supplier, chn, molclass, platebarcode, samplebarcode, storage)
            mollistsmiles = fast_fp_search(mollist, smiles, tanimoto)
            mollist = properties_search(mollistsmiles, minmw, maxmw, minlogp, maxlogp, minhba, maxhba, minhbd, maxhbd, mintpsa, maxtpsa, minfsp3, maxfsp3)
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
            mollist = id_search(cas, name, storageid, supplierid, supplier, chn, molclass, platebarcode, samplebarcode, storage)
            mollistsmart = smarts_search(mollist, smarts)
            mollist = properties_search(mollistsmart, minmw, maxmw, minlogp, maxlogp, minhba, maxhba, minhbd, maxhbd, mintpsa, maxtpsa, minfsp3, maxfsp3)
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
        mollist = id_search(cas, name, storageid, supplierid, supplier, chn, molclass, platebarcode, samplebarcode, storage)
        if cas == '':
            mollist = properties_search(mollist, minmw, maxmw, minlogp, maxlogp, minhba, maxhba, minhbd, maxhbd, mintpsa, maxtpsa, minfsp3, maxfsp3)
        
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
            #c += 1 # uncomment for testing
        except:
            pass
        #if c == 100: # uncomment for testing
        #    break # uncomment for testing
    logpmw = dumps(zip(logp, mw))
    hbamw = dumps(zip(hba, mw))
    hbdmw = dumps(zip(hbd, mw))
    tpsamw = dumps(zip(tpsa, mw))
    #histograms
    mwhist, mwedge = histogram(mw, bins=range(0,800,50))
    mwhist = dumps(zip(array(mwedge[:-1]).tolist(), array(mwhist).tolist()))
    fsp3hist, fsp3edge = histogram(fsp3, bins=20)
    fsp3hist = dumps(zip(array(fsp3edge[:-1]).tolist(), array(fsp3hist).tolist()))
    nrbhist, nrbedge = histogram(nrb, bins=[x + 0.5 for x in range(-1, 26, 1)])
    nrbhist = dumps(zip(array(nrbedge[:-1]).tolist(), array(nrbhist).tolist()))
    complexityhist, complexityedge = histogram(complexity, bins=range(0,160,10))
    complexityhist = dumps(zip(array(complexityedge[:-1]).tolist(), array(complexityhist).tolist()))
    logphist, logpedge = histogram(logp, bins=[x + 0.5 for x in range(-10, 10, 1)])
    logphist = dumps(zip(array(logpedge[:-1]).tolist(), array(logphist).tolist()))
    hbahist, hbaedge = histogram(hba, bins=[x + 0.5 for x in range(-1, 21, 1)])
    hbahist = dumps(zip(array(hbaedge[:-1]).tolist(), array(hbahist).tolist()))
    hbdhist, hbdedge = histogram(hbd, bins=[x + 0.5 for x in range(-1, 11, 1)])
    hbdhist = dumps(zip(array(hbdedge[:-1]).tolist(), array(hbdhist).tolist()))
    tpsahist, tpsaedge = histogram(tpsa, bins=range(0, 250, 10))
    tpsahist = dumps(zip(array(tpsaedge[:-1]).tolist(), array(tpsahist).tolist()))
    return render(request, 'statistics.html', {'servername':servername, 'logpmw':logpmw, 'hbamw':hbamw, 'hbdmw':hbdmw, 'tpsamw':tpsamw, 'mwhist':mwhist, 'fsp3hist':fsp3hist, 'nrbhist':nrbhist, 'complexityhist':complexityhist, 'logphist':logphist, 'hbahist':hbahist, 'hbdhist':hbdhist, 'tpsahist':tpsahist})

def upload(request):
    """
    Upload form with user handling. Can enable/disable functionality with variables at the top of this file.
    Set them to false to disable.
    """
    if 'confirm' in request.POST:
        mollist = Molecule.objects.filter(randomstring__contains=request.session.session_key)
        nummols = len(mollist)
        for mol in mollist:
            mol.randomstring = ""
            mol.save()
        success = []
        success.append(str(nummols) + " molecule(s) were added to the database!")
        return render(request, 'uploadresult.html', {'servername':servername, 'success':success})
    
    if 'cancel' in request.POST:
        Molecule.objects.filter(randomstring__contains=request.session.session_key).delete()
        success = []
        success.append("No molecules were added to the database!")
        return render(request, 'uploadresult.html', {'servername':servername, 'success':success})
    if 'logout' in request.POST:
        # Log out the user
        auth.logout(request)
        return HttpResponseRedirect(servername + "/upload.html")
    if 'single' in request.POST:
        form = SubmitSingle(request.POST)
        form.is_valid()
        error = []
        badsmiles = False
        randomstring = request.session.session_key #sessionid
        smiles = form.cleaned_data['smiles']
        name = form.cleaned_data['name']
        altname = form.cleaned_data['altname']
        cas = form.cleaned_data['cas']
        storage = form.cleaned_data['storage']
        storageID = form.cleaned_data['storageid']
        platebarcode = form.cleaned_data['platebarcode']
        samplebarcode = form.cleaned_data['samplebarcode']
        supplier = form.cleaned_data['supplier']
        supplierID = form.cleaned_data['supplierid']
        molclass = form.cleaned_data['molclass']
        comment = form.cleaned_data['comment']
        unit = form.cleaned_data['unit']
        amount = form.cleaned_data['amount']
        
        try:
            pybel.readstring("smi", str(smiles))
        except:
            error.append("Error in SMILES")
            badsmiles = True
        
        if form.is_valid() == False or badsmiles:
            for item in form.errors:
                error.append(form[item].errors)
            return render(request, 'uploadresult.html', {'servername':servername, 'error':error})
        
        addsingle(name, altname, supplier, supplierID, storage, storageID, unit, amount, cas, smiles, comment, molclass, platebarcode, samplebarcode, randomstring)
        mollist = Molecule.objects.filter(randomstring__contains=randomstring)
        return render(request, 'uploadresult.html', {'servername':servername, 'debugname':name, 'mollist':mollist})
    
    if 'table' in request.POST:
        form = UploadFileForm(request.POST, request.FILES)
        form.is_valid()
        error = []
        randomstring = request.session.session_key #sessionid
        if form.is_valid():
            cmdline = False
            userinput="yes"
            csvfile = request.FILES['file']#.read()
            csvdata = csv.reader(csvfile, dialect="excel")
            addtable(csvdata, cmdline, userinput, randomstring)
            mollist = Molecule.objects.filter(randomstring__contains=randomstring)
            return render(request, 'uploadresult.html', {'servername':servername, 'mollist':mollist})
        else:
            for item in form.errors:
                error.append(form[item].errors)
                return render(request, 'uploadresult.html', {'servername':servername, 'error':error})
    
    if uploadsingle == False and uploadtable == False:
        return render(request, 'upload.html', {'servername':servername, 'nogo':'nogo'})
    else:
        if logintoupload == True and not request.user.is_authenticated():
            # Serve notice that you need to log in in order to upload
            form = SubmitSingle()
            fileform = UploadFileForm()
            return render(request, 'upload.html', {'servername':servername, 'form':form, 'fileform':fileform})
        else:
            form = SubmitSingle()
            fileform = UploadFileForm()
            return render(request, 'upload.html', {'servername':servername, 'logedin':'ok', 'single':uploadsingle, 'table':uploadtable, 'form':form, 'fileform':fileform})
        
def login(request):
    """
    Log in/ log out functionality
    """
    if 'logout' in request.POST:
        auth.logout(request)
        return HttpResponseRedirect(servername + "/upload.html")
    if request.user.is_authenticated():
        return render(request, 'login.html', {'servername':servername, 'logout':'logout'})
    else:
        if 'login' in request.POST:
            username = request.POST['username']
            password = request.POST['password']
            user = auth.authenticate(username=username, password=password)
            if user is not None:
                if user.is_active:
                    auth.login(request, user)
                    return HttpResponseRedirect(servername + "/upload.html")
            # Redirect to a success page.
                else:
                    pass #TODO
            # Return a 'disabled account' error message
            else:
                pass #TODO
        else:
            return render(request, 'login.html', {'servername':servername})
        # Return an 'invalid login' error message.
    
