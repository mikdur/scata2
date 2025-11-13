from scata2.methods.models import ScataMethod
from django.db import models
from django.forms import ModelForm
from .cluster import do_cluster


class ScataDummyMethod(ScataMethod):
    dummy = models.CharField("Dummy setting", default="foo",
                             choices={"foo": "Foo",
                                      "bar": "Bar"})
    num_seqs = models.IntegerField("Number of sequences", editable=False,
                                   default=0)
    num_refs = models.IntegerField("Number of reference sequences", editable=False,
                                   default=0)
    mean_len = models.IntegerField("Mean sequence length", editable=False,
                                   default=0)
    num_clusters = models.IntegerField("Number of clusters (unique sequences)",
                                       editable=False, default=0)



    def cluster(self):
        do_cluster(self)

class ScataDummyMethodForm(ModelForm):

    class Meta:
        model = ScataDummyMethod
        fields = ["dummy"]
