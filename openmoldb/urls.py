from django.conf.urls import patterns, include, url
from openmoldb.views import home, molecule, browse, search, statistics, upload
# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    url(r'^$', home),
    url(r'^molecule/(\d+)/$', molecule),
    url(r'^browse', browse),
    url(r'^search', search),
    url(r'^stats', statistics),
    url(r'^upload', upload),
    # Examples:
    # url(r'^$', 'openmoldb.views.home', name='home'),
    # url(r'^openmoldb/', include('openmoldb.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
)
