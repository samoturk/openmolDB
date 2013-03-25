## Dependencies

OpenMolDB should run on any OS that can run Django 1.4.x.
It was tested on current version of Arch Linux and Arch Linux ARM (yes it runs on Raspberry Pi)

* Python 2.7.x or 2.6.x
* Django 1.4.x or 1.5.x
* OpenBabel 2.3.x with Python bindings (OpenBabel also needs cairo library)
* Numpy 1.7.x
* Any database supported by Django (with python bindings)
* Any web server that can run Django (with wsgi support)

## Installation

`git clone https://github.com/samoturk/openmolDB.github`

Then deploy it as you would any other Django project.

Set servername variable in openmoldb/views.py to the name of your server.
Have a look at settings.py and set your database engine, time zone, etc.

## Bulk addition of molecules

Save a table in your favourite spreadsheet program in csv format in excel
dialect. Recognizes colums named with (first row):
name, SMILES, altname, storage, storageID, CAS, supplier, supplierID, amount,
unit, comment, molclass, platebarcode, samplebarcode

Run: `python2 manage.py shell`
Import add2db module: `from add2db import add2db`
Finally add mols: `add2db("your_csv_table.csv")`
