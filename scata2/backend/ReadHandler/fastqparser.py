
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from .qualseq import Qual, QualSeq
from .exceptions import ScataReadsError, ScataFileError



class Single:
    
    def __init__(self, fastq_file):
        # Try to open as gzipped file, if that fails, open as plain fasta
        self.fastq = SeqIO.parse(fastq_file, "fastq")
        self.qual_present=True
        

    def __iter__ (self):
        return self

    def __next__(self):
        rec = next(self.fastq)
        q = Qual(rec.id, rec.letter_annotations["phred_quality"])
        rec.letter_annotations=dict()
        return QualSeq(rec, q)

class FastQPairQualSeq(QualSeq):
    def __init__(self, s1, s2, kmer, hsp, min):
        self.s1=s1
        self.s2=s2
        self.kmer = kmer
        self.hsp = hsp
        self.min = min
        self.s = None
        self.quals = None

    def get_seq(self):
        if not self.s:
            self._pair()
        return self.s
    
    def get_qual(self):
        if not self.s:
            self._pair()
        return self.quals
    

    def _pair(self):
        s1 = str(self.s1.seq)
        s2_rec = self.s2.reverse_complement()
        s2 = str(s2_rec.seq)

        # Build kmer table of read 2
        kmers2 = dict()
        for i in range(len(s2) - self.kmer):
            k = s2[i:i + self.kmer]
            kmers2.update({k : (kmers2.get(k, []) + [i])})
            
        # Identify read1 kmers in read2
        kmer_pos = [[x[0], kmers2.get(x[1]), x[2]] for x in
                     [(s1[y], s1[y : y + self.kmer], y) for y in range(len(s1) - self.kmer)]]

        # Identify the longest kmer run
        runs = [ ]
        r = [ ]
        for i in range(len(kmer_pos)):
            if kmer_pos[i][1]:
                r.append(kmer_pos[i])
            else:
                if len(r):
                    if len(r) > self.hsp:
                        runs.append(r)
                r=[]
        runs.sort(key=lambda a: len(a), reverse=True)
        

        if not len(runs) or sum([len(r) for r in runs]) < self.min:
            raise ScataReadsError("pairing_failed", "Pairing failed, no kmer runs / runs too short found")
        run = runs[0]

        # Remove ambigous start/end
        while len(run):
            if len(run[0][1]) > 1:
                run = run[1:]
            else:
                break
        while len(run):
            if len(run[-1][1]) > 1:
                run = run[:-1]
            else:
                break
        if not len(run):
            raise ScataReadsError("pairing_failed", "Pairing failed, no kmer run found")
        # Splice together new sequence
        try:
            new_seq = s1[:run[-1][2]] + s2[run[-1][1][0]:]
            quals = self.s1.letter_annotations["phred_quality"][:run[-1][2]] + \
                s2_rec.letter_annotations["phred_quality"][run[-1][1][0]:]
        except IndexError:
            print("IndexError")
            print(run, runs)
            raise ScataReadsError("pairing_failed", "Pairing failed, cpuld not splice")
        self.s=self.s2
        self.s.letter_annotations = dict()
        self.s.seq = Seq(new_seq)
        self.quals = Qual(self.s.id, quals)
                
class Pair:
    
    def __init__(self, fastq_1, fastq_2,
                 kmer=7, hsp=5, min=10):
        self.kmer=kmer
        self.hsp = hsp
        self.min = min
        try:
            handle1 = gzip.open(fastq_1, "rt")
            next(handle1) # This fails if not a gzip file, goes into the exception
            handle1 = gzip.open(fastq_1, "rt")
            self.fastq1 = SeqIO.parse(handle1, "fastq")
            self.fastq2 = SeqIO.parse(gzip.open(fastq_2, "rt"), "fastq")
        except IOError:
            self.fastq1 = SeqIO.parse(fastq_1, "fastq")
            self.fastq2 = SeqIO.parse(fastq_2, "fastq")
        self.qual_present = True


    def __iter__ (self):
        return self

    def __next__(self):
        s1 = next(self.fastq1)
        s2 = next(self.fastq2)

        return FastQPairQualSeq(s1, s2, self.kmer, self.hsp, self.min)
    

