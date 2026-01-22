from django.db import models
from django.urls import reverse
from django.contrib.auth.models import User
from django.core.validators import MinValueValidator, MaxValueValidator
from picklefield.fields import PickledObjectField
from .storages import get_file_storage, get_work_storage
from scata2.methods import methods as clustering_methods

import os.path

# Abstract base class for Scata models where the user interacts
# with the models


class ScataModel(models.Model):
    name = models.CharField("Name", max_length=50)
    description = models.TextField("Description", max_length=500,
                                   blank=True, null=True)
    errors = models.TextField("Error description", max_length=4000,
                              blank=True, default="")
    owner = models.ForeignKey(User, on_delete=models.CASCADE,
                              null=True, blank=True, editable=True)
    create_date = models.DateField("Create date", auto_now_add=True)
    deleted = models.BooleanField("Deleted", default=False)
    public = models.BooleanField("Public", default=False)

    def get_owner(self):
        d = ""
        if self.deleted:
            d = " (Deleted)"

        p = ""
        if self.public:
            p = " (Public)"

        if self.owner:
            return self.owner.username + d + p
        else:
            return "Orphaned" + d + p

    class Meta:
        abstract = True
        ordering = ['-create_date']


# Scata File class, used for uploaded files

class ScataFile(ScataModel):
    file = models.FileField("File", null=False,
                            storage=get_file_storage)
    file_size = models.PositiveBigIntegerField(default=0, editable=False)
    sha256 = models.CharField("Sha256 Sum", max_length=100, default="")

    def __str__(self):

        return "{u} {name} ({filename}) ({size}MB) {sha256}".\
            format(u=self.get_owner(),
                   name=self.name,
                   filename=os.path.basename(self.file.name),
                   size=self.file_size,
                   sha256=self.sha256)

    def get_absolute_url(self):
        return reverse("file-list")


# Scata Primer class. Used by ScataDataset

class ScataPrimer(ScataModel):
    sequence = models.CharField("Primer sequence (IUPAC ambiguity codes are " +
                                "accepted)", max_length=100)
    mismatches = models.DecimalField("Max mismatches",
                                     max_digits=2, decimal_places=0)

    def __str__(self):

        return "{u} {n}".format(u=self.get_owner(),
                                n=self.name)

    def get_absolute_url(self):
        return reverse("primer-list")


# Scata TagSet, used for demultiplexing
#

class ScataTagSet(ScataModel):
    is_valid = models.BooleanField("Valid", default=False,
                                   editable=False)
    validated = models.BooleanField("Validated", default=False,
                                    editable=False)
    num_tags = models.IntegerField("Number if tags", default=0,
                                   editable=False)
    tags = PickledObjectField(verbose_name="Data (dict)", editable=False,
                              null=True, compress=True)
    tagset_file = models.ForeignKey(ScataFile,
                                    null=True, on_delete=models.SET_NULL)

    def __str__(self):
        if self.is_valid and self.validated:
            status = ""
        elif not self.validated:
            status = " (Pending validation)"
        else:
            status = " (Failed, pending deletion)"

        if self.validated:
            size = " ({n} tags)".format(n=len(self.tags))
        else:
            size = ""

        return "{u} {name}{size}{status}".format(u=self.get_owner(),
                                                 name=self.name,
                                                 size=size,
                                                 status=status)

    def get_absolute_url(self):
        return reverse("tagset-list")

# Scata Amplicon (primers, tags and length limits)


class ScataAmplicon(ScataModel):
    five_prime_primer = models.ForeignKey(ScataPrimer, null=False, blank=False,
                                          on_delete=models.PROTECT,
                                          verbose_name="5' primer",
                                          related_name="five_prime_primer")
    five_prime_tag = models.ForeignKey(ScataTagSet, null=True,  blank=True,
                                       on_delete=models.PROTECT,
                                       verbose_name="5' tagset (not used " +
                                       "when clustering)",
                                       related_name="five_prime_tagset")
    three_prime_primer = models.ForeignKey(ScataPrimer, null=False,
                                           blank=False,
                                           on_delete=models.PROTECT,
                                           verbose_name="3' primer",
                                           related_name="three_prime_primer")
    three_prime_tag = \
        models.ForeignKey(ScataTagSet, null=True,  blank=True,
                          on_delete=models.PROTECT,
                          verbose_name="3' tagset (not used when clustering)",
                          related_name="three_prime_tagset")
    min_length = \
        models.IntegerField("Minimum length", default=0,
                            validators=[
                                MinValueValidator(0, "Minimum amplicon " +
                                                  "length is 0"),
                                MaxValueValidator(10000, "Max amplicon " + ""
                                                  "length is 10kbp")])
    max_length = \
        models.IntegerField("Maximum length (max 10kbp)",
                            default=10000,
                            validators=[MinValueValidator(0, "Minimum " +
                                                          "amplicon length " +
                                                          "is 0"),
                                        MaxValueValidator(10000, "Max " +
                                                          "amplicon length " +
                                                          "is 10kbp")])

    def __str__(self):
        t5 = self.five_prime_tag.name if self.five_prime_tag else ""
        p5 = self.five_prime_primer.name if self.five_prime_primer else ""
        p3 = self.three_prime_primer.name if self.three_prime_primer else ""
        t3 = self.three_prime_tag.name if self.three_prime_tag else ""

        return "{u} {name} ({t5}{p5} - {p3}{t3})".format(u=self.get_owner(),
                                                         name=self.name,
                                                         t5=t5, p5=p5,
                                                         p3=p3, t3=t3)

    def get_absolute_url(self):
        return reverse("amplicon-list")

# Set of reference sequences  with error model
#


class ScataReferenceSet(ScataModel):
    amplicon = \
        models.ForeignKey(ScataAmplicon, null=True, blank=True,
                          on_delete=models.PROTECT,
                          verbose_name="Filter database to contain region "
                          "defined by this amplicon. Leave empty to disable.")
    refseq_file = \
        models.ForeignKey(ScataFile, null=False, blank=False,
                          verbose_name="File with references",
                          on_delete=models.PROTECT)
    is_valid = models.BooleanField(default=False, editable=False)
    validated = models.BooleanField(default=False, editable=False)
    sequences = models.FileField("Sequences", null=True, blank=True,
                                 editable=False,
                                 upload_to="data/seqs",
                                 storage=get_work_storage)
    progress = models.CharField(default="", null=False, max_length=100)
    seq_count = models.IntegerField(editable=False, default=0)
    seq_total = models.IntegerField(editable=False, default=0)
    process_time = models.FloatField(editable=False, default=0.0)

    def __str__(self):
        if self.is_valid and self.validated:
            status = ""
        elif not self.validated:
            status = " (Pending validation)"
        else:
            status = " (Failed, pending deletion)"

        if self.validated:
            size = " ({n} sequences)".format(n=self.seq_count)
        else:
            size = ""

        return "{u} {name}{size}{status}".format(u=self.get_owner(),
                                                 name=self.name,
                                                 size=size,
                                                 status=status)

    def get_absolute_url(self):
        return reverse("referenceset-list")


class ScataRefsetErrorType(models.Model):
    refset = models.ForeignKey(ScataReferenceSet, on_delete=models.CASCADE)
    error = models.CharField(max_length=10)
    message = models.CharField(max_length=100)
    count = models.IntegerField()


# Dataset model, with error model and per tag statistics model
#

class ScataDataset(ScataModel):
    amplicon = models.ForeignKey(ScataAmplicon, null=True, blank=True,
                                 on_delete=models.PROTECT,
                                 verbose_name="Amplicon for filtering")
    mean_qual = \
        models.IntegerField("Minimum mean read quality", blank=False,
                            null=False, default=20,
                            validators=[MinValueValidator(10, "Min is 10"),
                                        MaxValueValidator(100, "Max is 100")])
    min_qual = \
        models.IntegerField("Minimum quality of any base", blank=False,
                            null=False, default=5,
                            validators=[MinValueValidator(1, "Min is 1"),
                                        MaxValueValidator(100, "Max is 100")])
    kmer_size = \
        models.IntegerField("Overlap join: kmer size", blank=False,
                            null=False, default=7,
                            validators=[MinValueValidator(5, "Min is 5"),
                                        MaxValueValidator(15, "Max is 15")])
    kmer_hsp_count = \
        models.IntegerField("Overlap join: HSP adjacent kmers", blank=False,
                            null=False, default=5,
                            validators=[MinValueValidator(3, "Min is 3"),
                                        MaxValueValidator(10, "Max is 10")])
    kmer_shared = \
        models.IntegerField("Overlap join: number of shared kmers",
                            blank=False, null=False, default=10,
                            validators=[MinValueValidator(5, "Min is 5"),
                                        MaxValueValidator(20, "Max is 20")])
    filter_method = \
        models.CharField("Filtering type", blank=False, null=False,
                         default="ampq", max_length=6,
                         choices={
                            "fsq": "Full sequence, quality screen",
                            "fs":  "Full sequence, NO quality filtering",
                            # "hqr":"Extract High Quality Region",
                            "ampq": "Amplicon quality",
                            "amp":  "Only amplicon extraction"})
    file_types = \
        models.CharField("File types", blank=False, null=False,
                         default="fastq", max_length=6,
                         choices={
                           "fastq":  "First file FASTQ",
                           "fastp":  "Paired End FASTQ (two files with mates)",
                           "fasta":  "First file FASTA",
                           "fastaq": "First file FASTA, second file quality " +
                                     "data"}
                         )

    # FIXME: Once the dataset is imported, the links to the files need
    # to be cleared. File name should be saved as string instead in order
    # to enable deletion of files in the future without loosing the dataset.

    file1 = models.ForeignKey(ScataFile, null=False, blank=False,
                              verbose_name="File 1",
                              on_delete=models.PROTECT,
                              related_name="dataset_file1")
    file2 = models.ForeignKey(ScataFile, null=True, blank=True,
                              verbose_name="File 2",
                              on_delete=models.PROTECT,
                              related_name="dataset_file2")

    is_valid = models.BooleanField(default=False, editable=False)
    validated = models.BooleanField(default=False, editable=True)
    has_stats = models.BooleanField(default=False, editable=False)

    progress = models.CharField(default="", null=False, max_length=100)
    seq_count = models.IntegerField(editable=False, default=0)
    seq_total = models.IntegerField(editable=False, default=0)
    seq_rev = models.IntegerField(editable=False, default=0)
    tag_count = models.IntegerField(editable=False, default=0)
    pc1_exp = models.FloatField(default=0)
    pc2_exp = models.FloatField(default=0)
    pc3_exp = models.FloatField(default=0)
    process_time = models.FloatField(editable=False, default=0.0)

    tags = models.FileField("Tag info", null=True, blank=True,
                            editable=False,
                            upload_to="data/tags",
                            storage=get_work_storage)

    sequences = models.FileField("Sequences", null=True, blank=True,
                                 editable=False,
                                 upload_to="data/seqs",
                                 storage=get_work_storage)

    def __str__(self):
        if self.is_valid and self.validated:
            status = ""
        elif not self.validated:
            status = " (Pending validation)"
        else:
            status = " (Failed, pending deletion)"

        if self.validated:
            size = " ({n} good sequences)".format(n=self.seq_count)
        else:
            size = ""

        return "{u} {name}{size}{status}".format(u=self.get_owner(),
                                                 name=self.name,
                                                 size=size,
                                                 status=status)

    def get_absolute_url(self):
        return reverse("dataset-list")


class ScataErrorType(models.Model):
    dataset = models.ForeignKey(ScataDataset, on_delete=models.CASCADE)
    error = models.CharField(max_length=100)
    message = models.CharField(max_length=100)
    count = models.IntegerField()


class ScataTagStat(models.Model):
    dataset = models.ForeignKey(ScataDataset, on_delete=models.CASCADE)
    tag = models.CharField(max_length=200)
    count = models.IntegerField()
    reversed = models.IntegerField()
    mean_len = models.FloatField()
    min_len = models.IntegerField()
    max_len = models.IntegerField()
    min_gc = models.FloatField()
    mean_gc = models.FloatField()
    max_gc = models.FloatField()
    in_pca = models.BooleanField(default=False)
    pc1 = models.FloatField(default=0)
    pc2 = models.FloatField(default=0)
    pc3 = models.FloatField(default=0)


# Scata job


class ScataJob(ScataModel):
    # Job settings
    amplicon = \
        models.ForeignKey("scata2.ScataAmplicon", null=True, blank=True,
                          on_delete=models.PROTECT,
                          verbose_name="Amplicon to use for clustering" +
                          "leave unset to use full imported amplicon")

    datasets = models.ManyToManyField(ScataDataset, verbose_name="Datasets",
                                      blank=False)
    refsets = models.ManyToManyField(ScataReferenceSet,
                                     verbose_name="References",
                                     blank=True)

    repseqs = models.IntegerField("Number of repseq to report", null=False,
                                  blank=False, default=3,
                                  validators=[MinValueValidator(0, "Min 0"),
                                              MaxValueValidator(100,
                                                                "Max 100")])

    method = models.CharField("Clustering method", blank=False,
                              default="scata", max_length=10,
                              choices={k: v['description'] for k, v in
                                       clustering_methods.items()})

    # Status and internal fields
    status = models.CharField("Status", default="Pending",
                              editable=False, max_length=250)

    def get_absolute_url(self):
        return reverse("job-list")

    def __str__(self):
        return "{u} {status}".format(u=self.get_owner(),
                                     status=self.status)

    def foo(self):
        print(clustering_methods)

# Global cluster
class ScataCluster(models.Model):
    job = models.ForeignKey(ScataJob, on_delete=models.CASCADE)

    size = models.IntegerField("Cluster size", null=False, blank=False, editable=False)
    num_reversed = models.IntegerField("Number of reversed sequences",
                                       null=False, blank=False, editable=False)
    num_genotypes = models.IntegerField("Number of genotypes", null=False, blank=False, editable=False)
    num_singletons = models.IntegerField("Number of singleton sequences",
                                         null=False, blank=False, editable=False)


class ScataClusterGenotype(models.Model):
    cluster = models.ForeignKey(ScataCluster, on_delete=models.CASCADE)
    sequence = models.TextField("Sequence", null=False, blank=False, editable=False)
    frequency = models.FloatField("Frequency", null=False, blank=False, editable=False)


class ScataTag(models.Model):
    job = models.ForeignKey(ScataJob, on_delete=models.CASCADE)
    size = models.IntegerField("Total reads", null=False, blank=False, editable=False)

class ScataTagCluster(models.Model):
    cluster = models.ForeignKey(ScataCluster, on_delete=models.CASCADE)
    tag = models.ForeignKey(ScataTag, on_delete=models.CASCADE)
    size = models.IntegerField("Cluster size", null=False, blank=False, editable=False)
    sequences = models.FileField("Cluster sequences", null=True, blank=True)
