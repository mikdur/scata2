from django.db import models
from django.forms import ModelForm


class DummyMethod(models.Model):
    job = models.ForeignKey("scata2.ScataJob", on_delete=models.CASCADE)
    dummy = models.CharField("Dummy setting", default="foo",
                             choices={"foo": "Foo",
                                      "bar": "Bar"})


class DummyMethodForm(ModelForm):

    class Meta:
        model = DummyMethod
        fields = ["dummy"]
