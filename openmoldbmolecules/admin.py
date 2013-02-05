from django.contrib import admin
from openmoldbmolecules.models import Molecule

class MoleculeAdmin(admin.ModelAdmin):
    list_display = ('id', 'name', 'added')

admin.site.register(Molecule, MoleculeAdmin)
