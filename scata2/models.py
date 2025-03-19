from django.db import models
from django.urls import reverse, reverse_lazy
from django.contrib.auth.models import User
from django.core.validators import MinValueValidator, MaxValueValidator
import json, os.path



# Abstract base class for Scata models where the user interacts
# with the models

class ScataModel(models.Model):
    owner = models.ForeignKey(User, on_delete = models.CASCADE,
                              null=True, blank=True, editable=True)
    create_date = models.DateField("Create date", auto_now_add=True)
    deleted = models.BooleanField("Deleted", default=False)
    public = models.BooleanField("Public", default=False)

    def get_owner(self):
        d=""
        if self.deleted:
            d = " (Deleted)"

        p=""
        if self.public:
            p=" (Public)"

        if self.owner:
            return self.owner.username + d + p
        else:
            return "Orphaned" + d + p

    
    class Meta:
        abstract = True


# Scata File class, used for uploaded files 

class ScataFile(ScataModel):
    name = models.CharField("File Name", max_length=50)
    description = models.TextField("Description", default="",
                                   max_length=500)
    file = models.FileField("File", upload_to="files/", null=False)
    file_size = models.PositiveBigIntegerField(default=0, editable=False)
    sha256 = models.CharField("Sha256 Sum", max_length=100, default="")
    
    def __str__(self):

        return "{u} {name} ({size}MB) {sha256}".format(u=self.get_owner(),
                                          name = os.path.basename(self.file.name),
                                          size=self.file_size,
                                          sha256=self.sha256)
    
    def get_absolute_url(self):
        return reverse("file-list")
    
# Scata Primer class. Used by ScataDataset 

class ScataPrimer(ScataModel):
    name = models.CharField("Short name", max_length=50)
    sequence = models.CharField("Primer sequence", max_length=100)
    mismatches = models.DecimalField("Max missmatches", 
                                     max_digits=2, decimal_places=0)
    description = models.TextField(max_length="1000")



    def __str__(self):
                
        return "{u} {n}".format(u=self.get_owner(), 
                                    n=self.name)
    
    def get_absolute_url(self):
        return reverse("primers-list")


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
    tagset_file = models.ForeignKey(ScataFile, 
                                    null=True, on_delete=models.SET_NULL)
    
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
    

class ScataAmplicon(ScataModel):
    name = models.CharField("Amplicon name", max_length=50)
    description = models.TextField("Description", max_length=300,
                                   blank=True, null=True)
    five_prime_primer = models.ForeignKey(ScataPrimer, null=True, blank=True,
                                          on_delete=models.SET_NULL, 
                                          verbose_name="5' primer",
                                          related_name="five_prime_primer")
    five_prime_tag = models.ForeignKey(ScataTagSet, null=True,  blank=True,
                                       on_delete=models.SET_NULL,
                                       verbose_name="5' tagset",
                                       related_name="five_prime_tagset")
    three_prime_primer = models.ForeignKey(ScataPrimer, null=True,  blank=True,
                                           on_delete=models.SET_NULL,
                                           verbose_name="3' primer",
                                           related_name="three_prime_primer")
    three_prime_tag = models.ForeignKey(ScataTagSet, null=True,  blank=True,
                                        on_delete=models.SET_NULL,
                                        verbose_name="3' tagset",
                                        related_name="three_prime_tagset")
    min_length=models.IntegerField("Minimum length", default=0,
                                   validators=[MinValueValidator(0, "Minimum amplicon length is 0"),
                                               MaxValueValidator(10000, "Max amplicon length is 10kbp")])
    max_length=models.IntegerField("Maximum length (max 10kbp)", default=10000,
                                   validators=[MinValueValidator(0, "Minimum amplicon length is 0"),
                                               MaxValueValidator(10000, "Max amplicon length is 10kbp")])
    
    def __str__(self):

        t5=self.five_prime_tag.name if self.five_prime_tag else ""
        p5=self.five_prime_primer.name if self.five_prime_primer else ""
        p3=self.three_prime_primer.name if self.three_prime_primer else ""
        t3=self.three_prime_tag.name if self.three_prime_tag else ""

        return "{u} {name} ({t5}{p5} - {p3}{t3})".format(u=self.get_owner(),
                                                   name = self.name,
                                                   t5=t5, p5=p5,
                                                   p3=p3, t3=t3)

    def get_absolute_url(self):
        return reverse("amplicon-list")



