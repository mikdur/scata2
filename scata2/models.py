from django.db import models
from django.contrib.auth.models import User


# Abstract base class for Scata models

class ScataModel(models.Model):
    owner = models.ForeignKey(User, on_delete = models.CASCADE,
                              null=True, blank=True, editable=True)
    create_date = models.DateField("Create date")
    public = models.BooleanField("Public", default=False)

    class Meta:
        abstract = True

class ScataPrimer(ScataModel):
    short_name = models.CharField("Short name", max_length=20)
    sequence = models.CharField("Primer sequence", max_length=100)
    missmatches = models.DecimalField("Max missmatches", 
                                     max_digits=2, decimal_places=0)
    description = models.TextField(max_length="1000")

    def __str__(self):
        if self.owner:
            user = self.owner.username
        else:
            user = "Orphaned"
        
        if self.public:
            public = "(Public)"
        else:
            public = ""

        return "{u} {p} {n}".format(u=user, p=public, n=self.short_name)

