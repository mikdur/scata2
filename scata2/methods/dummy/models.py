from scata2.methods.models import Method
from django.db import models
from django.forms import ModelForm
from .cluster import do_cluster


class DummyMethod(Method):
    dummy = models.CharField("Dummy setting", default="foo",
                             choices={"foo": "Foo",
                                      "bar": "Bar"})

    def cluster(self):
        do_cluster(self)
        print("done")
        pass


class DummyMethodForm(ModelForm):

    class Meta:
        model = DummyMethod
        fields = ["dummy"]
