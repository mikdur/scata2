# Parser to read and quality-filter Read results

from Bio import SeqIO
from Bio.Seq import Seq
from .qualseq import QualSeq, QualFile
from .filterseq import filter_full, filter_hqr, SeqDeTagger
from .fastqparser import Single, Pair
from .exceptions import ScataFileError, ScataReadsError




class RawReads:
    def __init__(self, fasta, qual_file=None):
        self.qual_present = True
        if qual_file:
            try:
                self.qual = QualFile(qual_file)
            except IOError:
                ScataFileError("bad_qualfile", "Bad .qual file")
        else:
            self.qual_present = False
        self.fasta = SeqIO.parse(fasta, "fasta")

    def __iter__ (self):
        return self

    def __next__(self):
        seq_record = next(self.fasta)
        qual = None
        if self.qual_present:
            qual = next(self.qual)
            if seq_record.id != qual.name:
                raise ScataFileError("id_mismatch", "Fasta and Qual ID missmatch")

            if len(seq_record.seq) != len(qual.quals):
                raise ScataFileError("length_mismatch", "Length of sequence and quality mismatch: " + seq_record.id + " " + qual.name)
        return QualSeq(seq_record, qual)






class Reads:
    def __init__(self, file1, file2=None,
                 file_type = "fastq", 
                 mean_min=20, min_qual=20, 
                 filtering="ampq",
                 amplicon=None,
                 kmer=7, hsp=5, hsp_min=10,
                 keep_primer=True,
                 ignore_tags = False):
        self.mean_min = int(mean_min)
        self.min_qual = int(min_qual)
        self.stats = dict(count = 0,
                          skipped = 0,
                          too_short = 0,
                          low_mean = 0,
                          low_min_quality = 0)
        self.filtering = filtering

        if amplicon:
            self.min_length = amplicon.min_length
            self.max_length = amplicon.max_length

            if ignore_tags:
                self.detagger = SeqDeTagger(amplicon,
                                        None,
                                        None,
                                        keep_primer=keep_primer)
            else:
                self.detagger = SeqDeTagger(amplicon,
                                            amplicon.five_prime_tag.tags if amplicon.five_prime_tag else None,
                                            amplicon.three_prime_tag.tags if amplicon.three_prime_tag else None,
                                            keep_primer=keep_primer)
    
        else:
            self.detagger = None
            self.min_length = 0
            self.max_length = 100000 # Arbitrary....



        if file_type == "fastq":
            self.rawreads = Single(file1)
        elif file_type == "fastq":
            self.rawreads = Pair(file1, file2,
                                 kmer=kmer, hsp=hsp, min=hsp_min)
        elif file_type == "fasta":
            self.rawreads = RawReads(file1)
        else:
            self.rawreads = RawReads(file1, file2)

    def __iter__ (self):
        return self

    def __next__ (self):

        qualseq = next(self.rawreads)
        
        if self.filtering == "fs":
            if self.detagger:
                return self.detagger.detag_seq(qualseq)
            else:
                return qualseq

        if not qualseq.qual:
            raise ScataFileError("missing_qual", "Selected filtering method requires quality data")
        
        qs = None
        if self.filtering == "fsq":
            qs = filter_full(qualseq, self.min_length,
                               self.mean_min, self.min_qual)
            if self.detagger:
                return self.detagger.detag_seq(qs)
            else:
                return qs

        elif self.filtering == "hqr":
            qs = filter_hqr(qualseq, self.min_length, 
                                  self.mean_min, self.min_qual)
            if self.detagger:
                return self.detagger.detag_seq(qs)
            else:
                return qs
        elif self.filtering == "ampq":
            if not self.detagger:
                raise ScataFileError("no_amplicon", "Amplicon quality reuires a defined amplicon")
            detagged = self.detagger.detag_seq(qualseq)
            return filter_full(detagged, self.min_length,
                               self.mean_min, self.min_qual)        


