from django import forms
import pybel

class SearchForm(forms.Form):
    moltext = forms.CharField(required=False, widget=forms.Textarea(attrs={'class':'span4', 'placeholder':'Mol data from editor', 'rows':'1'}))
    smiles = forms.CharField(required=False, widget=forms.TextInput(attrs={'class':'span4', 'placeholder':'SMILES'}))
    tanimoto = forms.FloatField(required=False, max_value=1.0, min_value=0.0, initial=0.8, error_messages={'max_value': 'Tanimoto: enter a number between 0 and 1!', 'min_value': 'Tanimoto: enter a number between 0 and 1!', 'invalid': 'Tanimoto: enter a number between 0 and 1!'}, widget=forms.TextInput(attrs={'class':'span1'}))
    minmw = forms.FloatField(required=False, initial=0, error_messages={'invalid': 'Molecular weight: enter a number!'}, widget=forms.TextInput(attrs={'class':'span1'}))
    maxmw = forms.FloatField(required=False, initial=9999, error_messages={'invalid': 'Molecular weight: enter a number!'}, widget=forms.TextInput(attrs={'class':'span1'}))
    minlogp = forms.FloatField(required=False, initial=-99, error_messages={'invalid': 'LogP: enter a number!'}, widget=forms.TextInput(attrs={'class':'span1'}))
    maxlogp = forms.FloatField(required=False, initial=99, error_messages={'invalid': 'LogP: enter a number!'}, widget=forms.TextInput(attrs={'class':'span1'}))
    minhba = forms.IntegerField(required=False, initial=0, error_messages={'invalid': 'H-bond acceptors: enter a whole number!'}, widget=forms.TextInput(attrs={'class':'span1'}))
    maxhba = forms.IntegerField(required=False, initial=99, error_messages={'invalid': 'H-bond acceptors: enter a whole number!'}, widget=forms.TextInput(attrs={'class':'span1'}))
    minhbd = forms.IntegerField(required=False, initial=0, error_messages={'invalid': 'H-bond donors: enter a whole number!'}, widget=forms.TextInput(attrs={'class':'span1'}))
    maxhbd = forms.IntegerField(required=False, initial=99, error_messages={'invalid': 'H-bond donors: enter a whole number!'}, widget=forms.TextInput(attrs={'class':'span1'}))
    chn = forms.CharField(required=False, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'CHNO'}))
    name = forms.CharField(required=False, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'Benzene'}))
    storageid = forms.CharField(required=False, max_length=10, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'xxxxx'}))
    cas = forms.CharField(required=False, max_length=12, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'xxxx-xx-x'}))
    supplier = forms.CharField(required=False, max_length=20, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'xxxxx'}))
    supplierid = forms.CharField(required=False, max_length=20, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'xxxxx'}))
    molclass = forms.CharField(required=False, max_length=20, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'chemical'}))

class SubmitSingle(forms.Form):
    name = forms.CharField(required=True, widget=forms.TextInput(attrs={'class':'span4', 'placeholder':'Benzene', 'required':''}))
    smiles = forms.CharField(required=True, widget=forms.TextInput(attrs={'class':'span4', 'placeholder':'SMILES', 'required':''}))
    altname = forms.CharField(required=False, widget=forms.TextInput(attrs={'class':'span4', 'placeholder':'Alt. name'}))
    storageid = forms.CharField(required=False, max_length=10, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'xxxxx'}))
    cas = forms.CharField(required=False, max_length=12, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'xxxx-xx-x'}))
    supplier = forms.CharField(required=False, max_length=20, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'xxxxx'}))
    supplierid = forms.CharField(required=False, max_length=20, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'xxxxx'}))
    molclass = forms.CharField(required=False, max_length=20, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'chemical'}))
    unit = forms.CharField(required=False, max_length=20, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'unit'}))
    amount = forms.CharField(required=False, max_length=20, widget=forms.TextInput(attrs={'class':'span2', 'placeholder':'xx'}))
    comment = forms.CharField(required=False, max_length=20, widget=forms.Textarea(attrs={'class':'span4', 'rows':'3', 'placeholder':'Comment'}))

class UploadFileForm(forms.Form):
    fileform  = forms.FileField(required=True)