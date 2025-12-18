import random
from time import sleep
from django.db import models


from scata2.methods.models import ScataMethod, ScataSequenceChunk, ChunkFullException
from django.forms import ModelForm
from django.core.validators import MinValueValidator, MaxValueValidator
import django_q.tasks as q2



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

        # Don't duplicate chunk set if already saved.

        chunks = list(ScataSequenceChunk.objects.filter(job=self.job).order_by("-length"))
        tasks = []

        if len(chunks) == 0:
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

        self.job.status = "Starting clustering"
        self.job.save()

        # In cases whith non-zero gap-penalty, clustering can be optimised by
        # only clustering sequences length difference less than
        # gap penalty * sequence length of the longer sequence. We use
        # global alignment, so length difference will always require a gap.

        # TODO implement the above. This is a plain all against all
        # implementation.

        task_group = "scata_cluster_{}".format(self.job.pk)
        tasks = []

        # Collect set of chunk groups of similar sequence counts
        chunks = list(ScataSequenceChunk.objects.filter(job=self.job).order_by("-length"))
        group_size = 100
        chunk_groups = []
        group = []
        c = 0
        while True:
            try:
                chunk = chunks.pop(0)
                c += len(chunk)
                group.append(chunk.pk)
            except IndexError:
                chunk_groups.append(group)
                break
            if c > group_size:
                chunk_groups.append(group)
                group = []
                c = 0

        print(chunk_groups)
        for q in range(len(chunk_groups)):
            for t in range(q, len(chunk_groups)):
                tasks.append(q2.async_task(ScataScataMethod.cluster_chunk,
                                           chunk_groups[q], chunk_groups[t],
                                           group=task_group,
                                           task_name="cluster_chunk self job={} {} {}". \
                                           format(self.job.pk,
                                                  len(chunk_groups[q]),
                                                  len(chunk_groups[t]))
                                           ))

        # Count groups while waiting.
        while True:
            success_count = q2.count_group(task_group)
            fail_count = q2.count_group(task_group, failures=True)
            total_count = success_count + fail_count

            self.job.status = "Clustering {}/{}".format(total_count,
                                              len(tasks))
            self.job.save()


            if total_count == len(tasks):
                break
            sleep(2)

        self.job.status = "Clustering done."
        self.job.save()

    @classmethod
    def cluster_chunk(cls, query, target):
        query = ScataSequenceChunk.objects.in_bulk(query)
        target = ScataSequenceChunk.objects.in_bulk(target)





class ScataScataMethodForm(ModelForm):

    class Meta:
        model = ScataScataMethod
        fields = ["distance", "min_alignment", "mismatch_pen",
                  "open_pen", "extend_pen", "endgap_pen", "max_homopolymer",
                  "downsample", "lowfreq"]


