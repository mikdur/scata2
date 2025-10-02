from django.db import models
from django.forms import ModelForm
from django.core.validators import MinValueValidator, MaxValueValidator


class ScataMethod(models.Model):
    job = models.ForeignKey("scata2.ScataJob", on_delete=models.CASCADE,
                            unique=True)
    distance = models.FloatField("Clustering distance 0.001 < x < 0.10",
                                 null=False, blank=False, default=0.015,
                                 validators=[MinValueValidator(0.001,
                                                               "Min 0.001"),
                                             MaxValueValidator(0.10,
                                                               "Max 0.10")])
    min_alignment = models.FloatField("Minimum alignment 0.5 < x < 1.0. Set to"
                                      " 1.0 to force global clustering (if all"
                                      " references are extracted amplicons)",
                                      null=False,
                                      blank=False, default=0.9,
                                      validators=[MinValueValidator(0.5,
                                                                    "Min 0.5"),
                                                  MaxValueValidator(1.0,
                                                                    "Max 1")])
    mismatch_pen = models.FloatField("Mismatch penalty 0 < x < 100", 
                                     null=False, blank=False, default=1,
                                     validators=[MinValueValidator(0, "Min 0"),
                                                 MaxValueValidator(100, 
                                                                   "Max 100")])
    open_pen = models.FloatField("Gap open penalty 0 < x < 100", null=False,
                                 blank=False, default=1,
                                 validators=[MinValueValidator(0, "Min 0"),
                                             MaxValueValidator(100, "Max 100")
                                             ])
    extend_pen = models.FloatField("Gap extension penalty 0 < x < 100",
                                   null=False, blank=False, default=0,
                                   validators=[MinValueValidator(0, "Min 0"),
                                               MaxValueValidator(100, 
                                                                 "Max 100")])
    endgap_pen = models.FloatField("End gap weight 0 < x < 100", null=False,
                                   blank=False, default=0,
                                   validators=[MinValueValidator(0, "Min 0"),
                                               MaxValueValidator(100, 
                                                                 "Max 100")])
    max_homopolymer = models.FloatField("Collaps homopolymers longer than " +
                                        "this. 0 to disable", null=False,
                                        blank=False, default=3,
                                        validators=[MinValueValidator(0, 
                                                                      "Min 0"),
                                                    MaxValueValidator(10, 
                                                                      "Max 10")
                                                    ])
    downsample = models.FloatField("Downsample samples with more than this "
                                   "number of reads to this number. 0 to " +
                                   "disable.", null=False,
                                   blank=False, default=0,
                                   validators=[MinValueValidator(0, "Min 0"),
                                               MaxValueValidator(1e5, 
                                                                 "Max 1e5")])
    lowfreq = models.FloatField("Remove global low frequency genotypes "
                                "occuring less than this number of times.", 
                                null=False, blank=False, default=0,
                                validators=[MinValueValidator(0, "Min 0"),
                                            MaxValueValidator(100, "Max 100")])


class ScataMethodForm(ModelForm):

    class Meta:
        model = ScataMethod
        fields = ["distance", "min_alignment", "mismatch_pen", "open_pen",
                  "extend_pen", "endgap_pen", "max_homopolymer",
                  "downsample", "lowfreq" ]
