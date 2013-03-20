Quick and dirty installation on Arch Linux
============
## Dependencies
`pacman -Sy git mod_wsgi2 apache python2-django openbabel python2-openbabel python2-psycopg2 postgresql python2-numpy python2-cairo cairo`

## Add user for OpenMolDB
`useradd -m -g users -s /bin/bash openmoldb`
`passwd openmoldb` 

## Login as "openmoldb" user
`git clone https://github.com/samoturk/openmolDB.github`
Move OpenMolDB files in /home/openmoldb/public_html

Change permissions: `chmod -R a+xr /home/openmoldb`
Open ~/public_html/openmoldb/views.py and set servername variable according to your configuration
Open ~/public_html/openmoldb/settings.py and set database varibles according to your configuration
Run `python2 ~/public_html/manage.py syncdb`

## Apache configuration
Add to the end of /etc/httpd/conf/httpd.conf:
`LoadModule wsgi_module modules/mod_wsgi.so
Alias /static/ /home/openmoldb/public_html/static/

<Directory /home/openmoldb/public_html/static>
Order deny,allow
Allow from all
</Directory>

WSGIScriptAlias / /home/openmoldb/public_html/openmoldb/wsgi.py
WSGIPythonPath /home/openmoldb/public_html

<Directory /home/openmoldb/public_html/openmoldb>
<Files wsgi.py>
Order deny,allow
Allow from all
</Files>
</Directory>

StartServers 1
MinSpareServers 1
MaxSpareServers 3
MaxClients 50`

Start apache: `systemctl start httpd`

## Postgresql configuration
Will add later