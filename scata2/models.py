from django.db import models
from django.urls import reverse, reverse_lazy
from django.contrib.auth.models import User
import json



# Abstract base class for Scata models where the user interacts
# with the models

class ScataModel(models.Model):
    owner = models.ForeignKey(User, on_delete = models.CASCADE,
                              null=True, blank=True, editable=True)
    create_date = models.DateField("Create date", auto_now_add=True)
    deleted = models.BooleanField("Deleted", default=False)
    public = models.BooleanField("Public", default=False)

    def get_owner(self):
        if self.deleted:
            d = " (Deleted)"
        else:
            d = ""

        if self.owner:
            return self.owner.username + d
        else:
            return "Orphaned" + d

    
    class Meta:
        abstract = True

# Scata Primer class. Used by ScataDataset 

class ScataPrimer(ScataModel):
    short_name = models.CharField("Short name", max_length=20)
    sequence = models.CharField("Primer sequence", max_length=100)
    mismatches = models.DecimalField("Max missmatches", 
                                     max_digits=2, decimal_places=0)
    description = models.TextField(max_length="1000")



    def __str__(self):
                
        if self.public:
            public = "(Public)"
        else:
            public = ""

        return "{u} {p} {n}".format(u=self.get_owner(), 
                                    p=public, n=self.short_name)
    
    def get_absolute_url(self):
        return reverse("primer-list")


# Scata TagSet, used for demultiplexing
# 

class ScataTagSet(ScataModel):
    name = models.CharField("Tagset Name", max_length=50)
    is_valid = models.BooleanField("Valid", default=False, 
                                   editable=False)
    validated = models.BooleanField("Validated", default=False, 
                                    editable=False)
    num_tags = models.IntegerField("Number if tags", default=0, 
                                   editable=False)
    tags = models.JSONField(encoder= json.JSONEncoder,
                            decoder= json.JSONDecoder,
                            verbose_name="Data (dict)", editable=False,
                            null=True)
    tagset_file = models.FileField("Tagset File", 
                                   upload_to="tagsets/",
                                   null=True)
    
    def __str__(self):
        if self.is_valid and self.validated:
            status = ""
        elif not self.validated:
            status = " (Pending validation)"
        else:
            status = " (Failed, pending deletion)"

        if self.validated:
            size = " ({n} tags)".format(n=len(self.tags))
        else:
            size = ""

        return "{u} {name}{size}{status}".format(u=self.get_owner(),
                                                   name = self.name,
                                               size=size,
                                               status=status)
    
    def get_absolute_url(self):
        return reverse("tagset-list")
    