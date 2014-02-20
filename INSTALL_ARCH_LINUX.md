Quick and dirty installation on Arch Linux
============
With nginx and gunicorn

## Dependencies
`pacman -Sy git nginx screen python2-django openbabel python2-openbabel python2-psycopg2 postgresql python2-numpy python2-cairo cairo`

Install [gunicorn](https://aur.archlinux.org/packages/gunicorn/) from AUR.

## Add user for OpenMolDB
`useradd -m -g users -s /bin/bash openmoldb`

`passwd openmoldb` 

## Login as "openmoldb" user

`mkdir /home/openmoldb/openmoldb`
`cd /home/openmoldb/openmoldb`
`git clone https://github.com/samoturk/openmolDB.github .`

In ~/openmoldb/openmoldb/views.py set variables under "Some settings" to either False or True.
Open ~/openmoldb/openmoldb/settings.py and set database varibles according to your configuration and add gunicorn to apps:
```python
INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.admin',
    'openmoldbmolecules',
    'gunicorn'
)
```

Run `python2 manage.py syncdb`
Run screen and in screen `python2 manage.py run_gunicorn`
Now gunicorn is serving the page.

## Nginx configuration
`mkdir /var/log/nginx`
`mkdir /var/log/nginx/logs`
Add to modify /etc/nginx/nginx.conf to look like:

```
user root;
worker_processes  1;

error_log  /var/log/nginx/logs/error.log;

events {
    worker_connections  1024;
}


http {
    include       mime.types;
    default_type  application/octet-stream;
    sendfile        on;
    keepalive_timeout  65;

    server {
        listen       80;
        server_name  localhost;

        access_log  /var/log/nginx/logs/host.access.log;
    # proxy for gunicorn
        location / {
            proxy_pass http://127.0.0.1:8000;
            proxy_set_header Host $host;
            proxy_set_header X-Real-IP $remote_addr;
            proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
            }
    # serving static files
    location /static/ {
            root /home/openmoldb/openmoldb;
        }

         }

    }
```

Start nginx: `sudo systemctl start nginx`

## Postgresql configuration
Will add later
