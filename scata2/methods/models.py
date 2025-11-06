from django.db import models
import gzip
import pickle
from Bio.SeqRecord import SeqRecord
from scata2.backend.ReadHandler.filterseq import SeqDeTagger
from scata2.backend.ReadHandler.qualseq import QualSeq
from scata2.backend.ReadHandler.exceptions import ScataReadsError



# Helper function to open dataset
def open_dataset(dataset):
    with dataset.sequences.file.open(mode="rb") as f:
        with gzip.open(f, mode="rb") as gz:
            return iter(pickle.load(gz).items())

class SeqIterator():
    detagger = None
    datasets = None
    current_dataset = None
    errors = dict()
    error_cnt = 0

    def __init__(self, datasets, amplicon=None):
        self.datasets = iter(datasets.all())
        if amplicon is not None:
            self.detagger = SeqDeTagger(amplicon)

    def __iter__(self):
        return self

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



class Method(models.Model):
    job = models.OneToOneField("scata2.ScataJob", on_delete=models.CASCADE,
                               related_name="job_method")

    seqs = None

    # def __init__(self, *args, **kwargs):
    #     super(Method, self).__init__(args, kwargs)
    #     self.seqs = SeqIterator(self.job.datasets, amplicon=self.job.amplicon)

    def get_seq_iterator(self):
        return SeqIterator(self.job.datasets, self.job.amplicon)

    def cluster(self):
        pass


    @classmethod
    def run_job(cls, job):
        instance = cls.objects.get(job=job)
        instance.cluster()
