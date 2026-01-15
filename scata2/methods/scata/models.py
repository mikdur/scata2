import gzip
import os
import pickle
import random
import subprocess
from io import BytesIO
from time import sleep

from django.core.files import File
from django.db import models
from django.conf import settings
from Bio import SeqIO

from scata2.storages import get_work_storage
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
        group_size = 4000
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
        task_num = 0
        for q in range(len(chunk_groups)):
            for t in range(q, len(chunk_groups)):
                task_num += 1
                tasks.append(q2.async_task(ScataScataMethod.cluster_chunk,
                                           self.job.pk, task_num,
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
    def cluster_chunk(cls, job_pk, task_num,
                      query, target):
        cls_instance = cls.objects.get(job=job_pk)

        target_file = os.path.join(settings.SCRATCH_DIR,
                                   "t_{}_{}.fasta".format(job_pk, task_num))
        query_file = os.path.join(settings.SCRATCH_DIR,
                                   "q_{}_{}.fasta".format(job_pk, task_num))
        query = ScataSequenceChunk.objects.in_bulk(query)
        target = ScataSequenceChunk.objects.in_bulk(target)

        query_records = [r for q in query.values() for r in q.get_uniseqs()]
        target_records = [r for t in target.values() for r in t.get_uniseqs()]

        SeqIO.write(target_records, target_file, "fasta")
        SeqIO.write(query_records, query_file, "fasta")

        vsearch_fields = { "query": str,
                           "target": str,
                           "id0": float,
                           "qilo": lambda a: int(a) - 1,
                           "qihi": lambda a: int(a) - 1,
                           "tilo": lambda a: int(a) - 1,
                           "tihi": lambda a: int(a) - 1,
                           "ql": int,
                           "tl": int,
                           "tcov": lambda a: float(a) / 100.0,
                           "qcov": lambda a: float(a) / 100.0,
                           "mism": int,
                           "opens": int,
                           "exts": int,
                           "pairs": int,
                           "pv": int,
                           "alnlen": int, }

        process = subprocess.Popen([settings.VSEARCH_COMMAND,
                          "--threads", "1",
                          "--usearch_global", query_file,
                          "--db", target_file,
                          "--userout", "-",
                          "--id", "0.95",
                          "--userfields", "+".join(vsearch_fields.keys()),
                          ], stdout=subprocess.PIPE, text=True)

        vsearch_result = process.communicate()[0]

        clusters = { }

        for line in vsearch_result.splitlines():
            hit = {a[0]: vsearch_fields[a[0]](a[1]) for a in zip(vsearch_fields.keys(), line.split("\t"))}

            # Ignore self
            if hit["query"] == hit["target"]:
                continue

            # Check alignment coverage
            if min(hit["tcov"], hit["qcov"]) < cls_instance.min_alignment:
                continue

            # Divergent sites
            distance = (hit["pairs"] - hit["pv"]) * cls_instance.mismatch_pen

            # Gaps
            distance += hit["opens"] * cls_instance.open_pen
            distance += hit["exts"] * cls_instance.extend_pen

            distance = distance / float(max(hit["qihi"] - hit["qilo"], hit["tihi"] - hit["tilo"]) + 1)

            # Check if within clusterin distance
            if distance > cls_instance.distance:
                continue

            # Case 1, both query and target in cluster, this joins
            # two already existing clusters.
            if hit["query"] in clusters and hit["target"] in clusters:
                cq = clusters[hit["query"]]
                ct = clusters[hit["target"]]

                if cq == ct:
                    continue # Same cluster, both already in

                # Merge clusters and set new cluster for all members
                nc = cq | ct
                for id in nc:
                    clusters[id] = nc

            elif hit["query"] in clusters:
                clusters[hit["query"]].add(hit["target"])
                clusters[hit["target"]] = clusters[hit["query"]]

            elif hit["target"] in clusters:
                clusters[hit["target"]].add(hit["query"])
                clusters[hit["query"]] = clusters[hit["target"]]

            else:
                clusters[hit["query"]] = { hit["target"], hit["query"] }
                clusters[hit["target"]] = clusters[hit["query"]]

        # Merge to unique list of clusters
        unique_clusters = []
        for cluster in clusters.values():
            if cluster not in unique_clusters:
                unique_clusters.append(cluster)

        if len(unique_clusters) > 0:
            ScataScataSubCluster.make_subcluster(unique_clusters, cls_instance.job)


class ScataScataSubCluster(models.Model):
    job = models.ForeignKey("scata2.ScataJob", on_delete=models.CASCADE)
    file = models.FileField(upload_to="scata/methods/scata/subcluster/", null=True, blank=True,
                            storage=get_work_storage)
    level = models.PositiveIntegerField(default=0)

    clusters = []

    @classmethod
    def make_subcluster(cls, clusters, job, level=0):
        cls_instance = cls()
        cls_instance.clusters = clusters
        cls_instance.job = job
        cls_instance.level = level
        cls_instance.save()

    def save(self, **kwargs):
        super().save(**kwargs) # Create db object to get pk
        with BytesIO() as seq_file:
            with gzip.open(seq_file, "wb") as gz:
                pickle.dump(self.clusters, gz)
            seq_file.seek(0)
            name = "j{}/l{}c{}".format(self.job.pk, self.level, self.pk)
            self.file = File(seq_file, name=name)
            super().save(**kwargs)


class ScataScataMethodForm(ModelForm):

    class Meta:
        model = ScataScataMethod
        fields = ["distance", "min_alignment", "mismatch_pen",
                  "open_pen", "extend_pen", "endgap_pen", "max_homopolymer",
                  "downsample", "lowfreq"]


