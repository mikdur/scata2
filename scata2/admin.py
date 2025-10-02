from django.contrib import admin

# Register your models here.

import scata2.models
import scata2.methods.scata.models
import scata2.methods.dummy.models


admin.site.register(scata2.models.ScataPrimer)
admin.site.register(scata2.models.ScataTagSet)
admin.site.register(scata2.models.ScataAmplicon)
admin.site.register(scata2.models.ScataFile)
admin.site.register(scata2.models.ScataReferenceSet)
admin.site.register(scata2.models.ScataJob)
admin.site.register(scata2.models.ScataDataset)
admin.site.register(scata2.methods.scata.models.ScataMethod)
admin.site.register(scata2.methods.dummy.models.DummyMethod)
