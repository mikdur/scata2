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

        # In cases whith non-zero gap-penalty, clustering can be optimised by
        # only clustering sequences length difference less than
        # gap penalty * sequence length of the longer sequence. We use
        # global alignment, so length difference will always require a gap.

        # TODO implement the above. This is a plain all against all
        # implementation.

        chunks = list(ScataSequenceChunk.objects.filter(job=self.job).order_by("-length"))

        def pop_chunks_n_seqs(_chunks, _n):
            chunks_n = []
            c=0
            while True:
                try:
                    chunks_n.append(_chunks.pop(0))
                    c += len(chunks_n[-1])
                except IndexError:
                    return chunks_n
                if c > _n:
                    return chunks_n

        chunk_size = 100
        while True:
            query = pop_chunks_n_seqs(chunks, chunk_size)
            if len(query) == 0:
                break

            #print("cluster_chunk self job={} {} {}".format(self.job.pk, len(query), len(query)))
            tasks.append(q2.async_task(ScataScataMethod.cluster_chunk,
                                       [q.pk for q in query],
                                       [q.pk for q in query],
                                       group="scata_cluster_{}".format(self.job.pk),
                                       task_name="cluster_chunk job={} {} {}".\
                                        format(self.job.pk, len(query), len(query))
                                       ))
            target = []
            c = 0
            if len(chunks) == 0:
                break
            for chunk in chunks:
                c += len(chunk)
                target.append(chunk)
                if c > chunk_size:
                    #print("cluster_chunk job={} {} {}".format(self.job.pk, len(query), len(target)))
                    tasks.append(q2.async_task(ScataScataMethod.cluster_chunk,
                                               [q.pk for q in query],
                                               [t.pk for t in target],
                                               group="scata_cluster_{}".format(self.job.pk),
                                               task_name="cluster_chunk job={} {} {}".\
                                                format(self.job.pk, len(query), len(target))
                                               ))
                    target = []
                    c=0
            # Cluster any leftovers
            if len(target) > 0:
                print("cluster_chunk job={} {} {}".format(self.job.pk, len(query), len(target)))
                tasks.append(q2.async_task(ScataScataMethod.cluster_chunk,
                                           [q.pk for q in query],
                                           [t.pk for t in target],
                                           group = "scata_cluster_{}".format(self.job.pk),
                                           task_name="cluster_chunk job={} {} {}".\
                                            format(self.job.pk, len(query), len(target))
                                           ))

        # Count groups while waiting.


    @classmethod
    def cluster_chunk(cls, query, target):
        print("running chunk {} {}".format(query, target))



class ScataScataMethodForm(ModelForm):

    class Meta:
        model = ScataScataMethod
        fields = ["distance", "min_alignment", "mismatch_pen",
                  "open_pen", "extend_pen", "endgap_pen", "max_homopolymer",
                  "downsample", "lowfreq"]


