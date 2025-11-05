from django.db import models
import gzip
import pickle


# Helper function to open dataset
def open_dataset(dataset):
    with dataset.sequences.file.open(mode="rb") as f:
        with gzip.open(f, mode="rb") as gz:
            return iter(pickle.load(gz).items())

class SeqIterator():
    # Init with queryset of sequences

    datasets = None
    current_dataset = None

    def __init__(self, datasets, amplicon=None):
        self.datasets = iter(datasets.all())

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_dataset is None:
            self.current_dataset = open_dataset(next(self.datasets))


        try:
            return next(self.current_dataset)
        except StopIteration:
            self.current_dataset = open_dataset(next(self.datasets))
            return next(self.current_dataset)


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
