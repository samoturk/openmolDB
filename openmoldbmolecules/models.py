from django.db import models

# Create your models here.
#class molclass(models.Model):
#    nameOfType = models.CharField(max_length=20)
#    def __unicode__(self):
#        return self.nameOfType

class Molecule(models.Model):
    """
    Molecule models
    Only name needed, date is added automatically
    """
    name = models.CharField(max_length=500, blank=True) # User defined
    altname = models.CharField(max_length=500, blank=True) # User defined
    storageID = models.CharField(max_length=30, blank=True) # User defined
    storage = models.CharField(max_length=30, blank=True) # User defined TODO show in molecule model
    CAS = models.CharField(max_length=12, blank=True) # User defined
    supplier = models.CharField(max_length=100, blank=True) # User defined
    supplierID = models.CharField(max_length=100, blank=True) # User defined
    altsupplier = models.CharField(max_length=100, blank=True) # User defined, not shown
    altsupplierID = models.CharField(max_length=100, blank=True) # User defined, not shown
    SMILES = models.CharField(max_length=500, blank=True) # User defined, basis for calculated properties
    RMW = models.FloatField(null=True, blank=True, verbose_name='MW, provided by supplier') # Not implemented/used
    CMW = models.FloatField(null=True, blank=True, verbose_name='MW, calculated') # Calculated from SMILES
    CHN = models.CharField(max_length=20, blank=True, verbose_name='Brutto formula') # Calculated from SMILES
    HBA = models.IntegerField(null=True, blank=True, verbose_name='H-bond acceptors') # Calculated from SMILES
    HBD = models.IntegerField(null=True, blank=True, verbose_name='H-bond donors') # Calculated from SMILES
    logP = models.FloatField(null=True, blank=True) # Calculated from SMILES
    tpsa = models.FloatField(null=True, blank=True, verbose_name='Polar surface area') # Calculated from SMILES
    fingerprint = models.TextField(blank=True, verbose_name='Fingerprint') # Calculated from SMILES
    nrb = models.IntegerField(null=True, blank=True, verbose_name='Number of rotatable bonds') # Calculated from SMILES
    fsp3 = models.FloatField(null=True, blank=True) # Calculated from SMILES
    complexity = models.IntegerField(null=True, blank=True, verbose_name='Molecular complexity') # Calculated from SMILES
    amount = models.FloatField(null=True, blank=True, verbose_name='Quantity') # User defined
    unit = models.CharField(max_length=10, blank=True, verbose_name='Unit') # User defined
    comment = models.TextField(blank=True, verbose_name='Comment') # User defined
    molfile = models.TextField(blank=True, verbose_name='MolFile') # Calculated from SMILES
    pains = models.TextField(max_length=100, verbose_name='Pan Assay Interference Compounds', blank=True) # Calculated from SMILES
    added = models.DateField(auto_now=True, auto_now_add=True) #automatically updates date and adds date upon creation
    molclass = models.CharField(max_length=30, verbose_name='Type of compound') # User defined (recommended: chemical, hts, compound; for hts and compound PAINS can be calculated
    randomstring = models.CharField(max_length=32, blank=True) # Uses session id
    platebarcode = models.CharField(max_length=64, verbose_name='Plate barcode', blank=True) # User defined
    samplebarcode = models.CharField(max_length=64, verbose_name='Sample barcode', blank=True) # User defined
    
    def __unicode__(self):
        return u'%s, %s' % (self.name, self.CAS)
    
    class Meta:
        ordering = ['name']

