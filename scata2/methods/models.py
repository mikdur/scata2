from io import BytesIO
from django.db import models
from django.core.files import File
import gzip
import pickle
from Bio.SeqRecord import SeqRecord
from scata2.storages import get_work_storage
from scata2.backend.ReadHandler.filterseq import SeqDeTagger
from scata2.backend.ReadHandler.qualseq import QualSeq
from scata2.backend.ReadHandler.exceptions import ScataReadsError



# Helper function to open dataset
def open_dataset(dataset):
    with dataset.sequences.file.open(mode="rb") as f:
        with gzip.open(f, mode="rb") as gz:
            return iter(pickle.load(gz).items())

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
    def new_chunk(cls, job, length, chunk_size=2000):
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



class ScataOTU(models.Model):
    job = models.ForeignKey("scata2.ScataJob", on_delete=models.CASCADE)
    name = models.CharField(max_length = 200, default="")
    size = models.IntegerField(default = 0)


class ScataRepSeq(models.Model):
    otu = models.ForeignKey(ScataOTU, on_delete=models.CASCADE)
    sequence = models.CharField(max_length = 2000, default = "")
    length = models.IntegerField(default = 0)
    frequency = models.IntegerField(default = 0)





