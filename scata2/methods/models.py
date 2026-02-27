import math
from io import BytesIO
from django.db import models
from django.core.files import File
import gzip
import pickle
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from django.http import Http404

import numpy as np
from scipy.special import gammaln
from sklearn.cluster import AgglomerativeClustering

from scata2.storages import get_work_storage
from scata2.backend.ReadHandler.filterseq import SeqDeTagger
from scata2.backend.ReadHandler.qualseq import QualSeq
from scata2.backend.ReadHandler.exceptions import ScataReadsError


# Helper function to open dataset
def open_dataset(dataset):
    with dataset.sequences.file.open(mode="rb") as f:
        with gzip.open(f, mode="rb") as gz:
            return iter(pickle.load(gz).items())

def open_tags(dataset):
    with dataset.tags.file.open(mode="rb") as f:
        with gzip.open(f, mode="rb") as gz:
            return pickle.load(gz)

class SeqIterator():
    total_cnt = 0
    error_cnt = 0
    current_dataset = None
    amplicon = None
    detagger = None
    datasets = None
    errors = dict()

    def __init__(self, datasets, amplicon=None):
        for dataset in datasets.all():
            self.total_cnt += dataset.seq_count
        self.datasets = iter(datasets.all())
        if amplicon is not None:
            self.detagger = SeqDeTagger(amplicon)

    def __iter__(self):
        return self

    def __len__(self):
        return self.total_cnt

    def __next__(self):
        if self.current_dataset is None:
            self.current_dataset = open_dataset(next(self.datasets))

        if self.detagger is None:
            try:
                return next(self.current_dataset)
            except StopIteration:
                self.current_dataset = open_dataset(next(self.datasets))
                return next(self.current_dataset)
        else:
            while True:
                try:
                    seq_tuple = next(self.current_dataset)
                except StopIteration:
                    self.current_dataset = open_dataset(next(self.datasets))
                    seq_tuple = next(self.current_dataset)

                try:
                    seq = self.detagger.detag_seq(QualSeq(SeqRecord(seq_tuple[1]))).seq_record.seq
                    return (seq_tuple[0], seq)
                except ScataReadsError as e:
                    self.error_cnt += 1
                    if e.error not in self.errors:
                        self.errors[e.error] = dict( cnt = 1,
                                                     msg = e.message)
                    else:
                        self.errors[e.error]['cnt'] += 1


class RefIterator():
    total_cnt = 0
    current_refset = None
    refsets = None


    def __init__(self, refsets):
        for refset in refsets.all():
            self.total_cnt += refset.seq_count
        self.refsets = iter(refsets.all())

    def __len__(self):
        return self.total_cnt

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_refset is None:
            self.current_refset = open_dataset(next(self.refsets))

        try:
            return next(self.current_refset)
        except StopIteration:
            self.current_refset = open_dataset(next(self.refsets))
            return next(self.current_refset)


class ScataMethod(models.Model):

    job = models.OneToOneField("scata2.ScataJob", on_delete=models.CASCADE)

    seqs = None

    def get_seq_iterator(self):
        return SeqIterator(self.job.datasets, self.job.amplicon)

    def get_ref_iterator(self):
        return RefIterator(self.job.refsets)

    def cluster(self):
        pass


    @classmethod
    def run_job(cls, job):
        instance = cls.objects.get(job=job)
        instance.cluster()

    # Called from view to generate data for visualisation
    def get_facet(self, facet):

        if facet == "clusters":
            return self.get_clusters()
        elif facet == "clustertable":
            return self.get_clustertable()
        elif facet == "clustertag_relative":
            return self.get_clustertag_relative()
        elif facet == "species_accumulation":
            return self.get_species_accumulation()
        else:
            raise Http404("No such facet")


    # Return list of clusters and their sizes
    def get_clusters(self):
        clusters = ScataCluster.objects.filter(job=self.job).order_by("-size")[:200]
        return [{"id": c.name,
                 "size": c.size,
                 "size_log2": math.log2(c.size),
                 "size_log10": math.log10(c.size),
                 "singletons": c.num_singletons / c.size,
                 "genotypes": c.num_genotypes } for c in clusters]

    def get_clustertable(self):
        clusters = ScataCluster.objects.filter(job=self.job).order_by("-size")[:40]

        return [ { "name": c.name,
                   "size": c.size,
                   "genotypes": c.num_genotypes,
                   "singletons": c.num_singletons,
                   #"Reversed": c.num_reversed,
                 } for c in clusters ]


    def get_clustertag_relative(self):
        clusters = [c.pk for c in ScataCluster.objects.filter(job=self.job).order_by("-size")[:60]]
        tags = ScataTag.objects.filter(job=self.job).order_by("name")


        pk2i = { v:i for i,v in enumerate(clusters) }

        X = np.zeros((len(tags), len(clusters)), dtype=np.float64)

        ret = []
        for i,t in enumerate(tags):
            tagclusters = ScataTagCluster.objects.filter(tag=t, cluster__in=clusters)
            tag_size = sum([tc.size for tc in tagclusters])
            for tc in tagclusters:
                ret.append({ "tag": t.name,
                     "cluster": tc.cluster.name,
                     "size": "{}".format(tc.size),
                     "rel_size": "{}".format((tc.size / tag_size)), })
                X[i, pk2i[tc.cluster.pk]] = math.log2(tc.size / tag_size)

        #order = leaves_list(linkage(X, optimal_ordering=True))

        model = AgglomerativeClustering(distance_threshold=0, n_clusters=None)
        model.fit(X)

        order = [int(x) for c in model.children_ for x in c if x < model.n_leaves_]


        tagname2order = {}
        for i in range(len(order)):
            tagname2order[tags[i].name] = order[i]

        for i in range(len(ret)):
            ret[i]['c_order'] = int(tagname2order[ret[i]['tag']])

        return ret


    def get_species_accumulation(self):

        # https://stackoverflow.com/questions/26938888/log-computations-in-python
        # n choose k can be implemented using gamma distribution

        def _combln(n, k):
            return gammaln(n + 1) - gammaln(n - k + 1)

        tags = ScataTag.objects.filter(job=self.job).order_by("name")
        ret = []
        for tag in tags:
            stc = ScataTagCluster.objects.filter(tag=tag)

            sizes = np.array([a.size for a in stc], dtype=np.float64)

            K = len(sizes)
            N = sum(sizes)

            curve = [ ]
            for n in np.linspace(1, N, num=50, dtype=np.float64):
                summation = np.float64(0.0)
                for s in sizes:
                    summation += np.exp(_combln(N - s, n) - _combln(N, n))
                curve.append((float(n), float(K - summation)))

            ret += [{ "tag": tag.name,
                         "x": a[0],
                         "y": a[1]} for a in curve]
        return ret


# Models to represent chunk of sequences
class ChunkFullException(Exception):
    def __init__(self, error="File full", message="Chunk full"):
        self.error = error
        self.message = message
        super().__init__()


# Needs dereplication
class ScataSequenceChunk(models.Model):
    job = models.ForeignKey("scata2.ScataJob", on_delete=models.CASCADE)
    length = models.IntegerField(default = 0)
    num_sequences = models.IntegerField(default = 0)
    num_uniques = models.IntegerField(default = 0)
    errors = dict()
    file = models.FileField(upload_to = "scata/methods/scata/chunk/", null=False,
                            storage=get_work_storage)

    def __init__(self, *args, **kwargs):
        self.sequences = {}
        super().__init__(*args, **kwargs)

    def __str__(self):
        return "ScataSequenceChunk(job={}, l={}, n={}, u={})".format(self.job, self.length, self.num_sequences, self.num_uniques)

    def __len__(self):
        return self.num_uniques

    @classmethod
    def new_chunk(cls, job, length, chunk_size=4000):
        chunk = cls()
        chunk.job = job
        chunk.length = length
        chunk.chunk_size = chunk_size
        chunk.num_sequences = 0
        return chunk

    def add_sequence(self, id, sequence):
        assert(len(sequence) == self.length)
        self.num_sequences += 1
        self.sequences[str(sequence)] = self.sequences.get(str(sequence), []) + [id]
        self.num_uniques = len(self.sequences)
        if len(self.sequences) > self.chunk_size:
            raise(ChunkFullException())


    def save(self, **kwargs):
        with BytesIO() as seq_file:
            with gzip.open(seq_file, "wb") as gz:
                pickle.dump(self.sequences, gz)
            seq_file.seek(0)
            name = "j{}/s{}c{}".format(self.job.pk, self.length, self.num_sequences)
            self.file = File(seq_file, name=name)
            super().save(**kwargs)

    def _load(self):
        if len(self.sequences) == 0:
            self.sequences = {}
            with self.file.open(mode="rb") as f:
                with gzip.open(f, mode="rb") as gz:
                    self.sequences = pickle.load(gz)

    def get_uniseqs(self):
        self._load()

        return [SeqRecord( seq = Seq(a[1]), id="{}_{}".format(self.pk, a[0]),
                           description="")
                for a in enumerate(self.sequences.keys())]

    def get_seq_by_id(self, item):
        self._load()
        return self.sequences[list(self.sequences.keys())[item]]


class ScataOTU(models.Model):
    job = models.ForeignKey("scata2.ScataJob", on_delete=models.CASCADE)
    name = models.CharField(max_length = 200, default="")
    size = models.IntegerField(default = 0)


class ScataRepSeq(models.Model):
    otu = models.ForeignKey(ScataOTU, on_delete=models.CASCADE)
    sequence = models.CharField(max_length = 2000, default = "")
    length = models.IntegerField(default = 0)
    frequency = models.IntegerField(default = 0)


# Global cluster
class ScataCluster(models.Model):
    job = models.ForeignKey("scata2.ScataJob", on_delete=models.CASCADE)

    name = models.CharField(max_length = 100, default = "")

    size = models.IntegerField("Cluster size", null=False, blank=False, editable=False,
                               default=0)
    num_reversed = models.IntegerField("Number of reversed sequences",
                                       null=False, blank=False, editable=False,
                                       default=0)
    num_genotypes = models.IntegerField("Number of genotypes", null=False, blank=False, editable=False,
                                        default=0)
    num_singletons = models.IntegerField("Number of singleton sequences",
                                         null=False, blank=False, editable=False,
                                         default=0)


class ScataClusterGenotype(models.Model):
    cluster = models.ForeignKey(ScataCluster, on_delete=models.CASCADE)
    sequence = models.TextField("Sequence", null=False, blank=False, editable=False)
    frequency = models.FloatField("Frequency", null=False, blank=False, editable=False)


class ScataTag(models.Model):
    job = models.ForeignKey("scata2.ScataJob", on_delete=models.CASCADE)
    size = models.IntegerField("Total reads", null=False, blank=False, editable=False,
                               default=0)
    num_clusters = models.IntegerField("Number of clusters", null=False, blank=False, editable=False,
                                       default=0)
    name = models.CharField("Name", max_length=200, null=False, blank=False, editable=False,
                            default="")

class ScataTagCluster(models.Model):
    cluster = models.ForeignKey(ScataCluster, on_delete=models.CASCADE)
    tag = models.ForeignKey(ScataTag, on_delete=models.CASCADE)
    size = models.IntegerField("Cluster size", null=False, blank=False, editable=False,
                               default=0)
    sequences = models.FileField("Cluster sequences", upload_to = "scata/methods/scata/tag/",
                                 storage=get_work_storage, null=True, blank=True)

    seq_list = []


    def save(self, **kwargs):
        super().save(**kwargs)
        with BytesIO() as seq_file:
            with gzip.open(seq_file, "wb") as gz:
                pickle.dump(self.sequences, gz)
            seq_file.seek(0)
            name = "j{}/c{}_t{}".format(self.cluster.job.pk, self.cluster.pk, self.pk)
            self.file = File(seq_file, name=name)
            super().save(**kwargs)



