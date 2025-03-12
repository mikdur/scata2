from django.contrib import admin

# Register your models here.

import scata2.models

admin.site.register(scata2.models.ScataPrimer)
admin.site.register(scata2.models.ScataTagSet)

