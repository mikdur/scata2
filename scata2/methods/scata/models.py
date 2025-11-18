from django.db import models
from scata2.methods.models import ScataMethod, ScataSequenceChunk, ChunkFullException
from django.forms import ModelForm
from django.core.validators import MinValueValidator, MaxValueValidator


class ScataScataMethod(ScataMethod):
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
    max_homopolymer = models.FloatField("Collapse homopolymers longer than " +
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
                                "occurring less than this number of times.",
                                null=False, blank=False, default=0,
                                validators=[MinValueValidator(0, "Min 0"),
                                            MaxValueValidator(100, "Max 100")])

    # Clustering method
    def cluster(self):
        print("SCATA Clustering {}".format(self))
        self.job.status = "Preparing"
        self.job.save()

        seq_iter = self.get_seq_iterator()



        id2name = {}
        seqs = {}
        n = 0

        self.job.status = "Deduplicating 0/{}".format(len(seq_iter))
        self.job.save()

        for seq in seq_iter:
            n += 1
            if n % 10000 == 0:
                self.job.status = "Deduplicating {}/{}".format(n, len(seq_iter))
                self.job.save()
                print("Deduplicating {}/{}".format(n, len(seq_iter)))
            id = "s{}".format(n)
            l = len(seq[1])
            id2name[id] = seq[0]
            chunk = seqs.get(l, ScataSequenceChunk.new_chunk(self.job, l))
            try:
                chunk.add_sequence(id, seq[1])
            except ChunkFullException:
                chunk.save()
                old_chunk = chunk
                chunk = None
                chunk = ScataSequenceChunk.new_chunk(self.job, l)
                assert(chunk != old_chunk)
                chunk.add_sequence(id, seq[1])

            seqs[l] = chunk

        for seq in seqs.values():
            seq.save()

        # In cases whith non-zero gap-penalty, clustering can be optimised by
        # only clustering sequences length difference less than
        # gap penalty * sequence length of the longer sequence. We use
        # global alignment, so length difference will always require a gap.


class ScataScataMethodForm(ModelForm):

    class Meta:
        model = ScataScataMethod
        fields = ["distance", "min_alignment", "mismatch_pen",
                  "open_pen", "extend_pen", "endgap_pen", "max_homopolymer",
                  "downsample", "lowfreq"]


