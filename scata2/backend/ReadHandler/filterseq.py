from Bio.Seq import Seq
from .exceptions import ScataReadsError
from .qualseq import QualSeq


def filter_full(qualseq, min_length, mean_min, min_qual):
    seq_record = qualseq.get_seq()
    qual = qualseq.get_qual()

    if len(seq_record.seq) < min_length:
        raise ScataReadsError("too_short", "Read too short")
        
    if float(sum(qual.quals)) / len(qual.quals) < mean_min:
        raise ScataReadsError("low_mean_quality", "Too low mean quality")
            
    seq_list = list(seq_record.seq)
            
    seq_list = list(map(lambda b, q: b if q >= min_qual else 'N',
                   seq_list, qual.quals))
    
    if seq_list.count('N'):
        raise ScataReadsError("low_min_quality", "Too low minimum quality")
    
    seq_record.seq = Seq("".join(seq_list))
    
    qualseq.seq_record = seq_record

    return qualseq

# This code is broken. Leaving it to be fixed
# if there is demand for it.

def filter_hqr(qualseq, min_length, mean_min, min_qual):
    seq_record = qualseq.get_seq()
    qual = qualseq.get_qual()
    if len(seq_record.seq) < min_length:
        raise ScataReadsError("too_short", "Read too short")

    # Get regions with quality above minimum
                    
    hqr = [ ]
                    
    t_start=-1
    t_end=-1

    for i in range(len(qual.quals)):
        if t_start == -1: # Outside region
            if qual.quals[i] >= min_qual:
                t_start = i
        else: # Inside region
            if qual.quals[i] < min_qual:
                t_end = i
                hqr.append(dict(start = t_start,
                                end = t_end,
                                mean_quality = 0,
                                length = 0))
                t_start = -1
                t_end = -1
    if t_start != -1 and t_end == -1:
        hqr.append(dict(start = t_start,
                        end = len(qual.quals) - 1,
                        mean_quality = 0,
                        length = 0))
    for r in hqr:
        if r["end"] - r["start"] < min_length:
            continue
        while r["end"] > r["start"] and \
                (sum(qual.quals[r["start"]:r["end"]]) / float(r["end"] - r["start"])) \
                < mean_min:
            if qual.quals[r["start"]] < qual.quals[r["end"]]:
                r["start"] += 1
            else:
                r["end"] -= 1
        if r["end"] > r["start"]:
            r["mean_qual"] = (sum(qual.quals[r["start"]:r["end"]]) / float(r["end"] - r["start"]))
            r["length"] = r["end"] - r["start"]

    hqr = [a for a in hqr if a["length"] >= min_length and a["mean_qual"] >= mean_min]
    
    if not len(hqr):
        raise ScataReadsError("low_mean_quality", "Too low mean quality")

    hqr.sort(key=lambda a: a["length"], reverse=True)

                    
    if hqr[0]["end"] - hqr[0]["start"] < min_length:
        raise ScataReadsError("low_mean_quality", "Too low mean quality")
                    
    seq_record.seq = seq_record.seq[hqr[0]["start"]:hqr[0]["end"]]
    qual.quals = qual.quals[hqr[0]["start"]:hqr[0]["end"]]

    qualseq.seq_record = seq_record
    qualseq.qual = qual
    return qualseq


class DeTaggedSeq(QualSeq):
    def __init__(self, seq_record=None, qual=None, tag="", reversed=False):
        self.seq_record = seq_record
        self.qual = qual
        self.tag = tag
        self.reversed = reversed

    def __str__(self):
        return "DeTaggedSeq({name}, tag={tag}, len={len})".format(
            name=self.seq_record.id,
            tag=self.tag if self.tag else "None",
            len=len(self.seq_record.seq)
        )        

base_A = 1
base_C = 2
base_G = 4
base_T = 8
trans_table = { 'A' : base_A,
                'C' : base_C,
                'G' : base_G,
                'T' : base_T,
                'R' : base_A | base_G,
                'Y' : base_C | base_T,
                'S' : base_G | base_C,
                'W' : base_A | base_T,
                'K' : base_G | base_T,
                'M' : base_A | base_C,
                'B' : base_C | base_G | base_T,
                'D' : base_A | base_G | base_T,
                'H' : base_A | base_C | base_T,
                'V' : base_A | base_C | base_G,
                'N' : base_A | base_C | base_T | base_G }

complement_trans_table = { 'A' : base_T,
                'C' : base_G,
                'G' : base_C,
                'T' : base_A,
                'R' : base_T | base_C,
                'Y' : base_G | base_A,
                'S' : base_C | base_G,
                'W' : base_T | base_A,
                'K' : base_C | base_A,
                'M' : base_T | base_G,
                'B' : base_G | base_C | base_A,
                'D' : base_T | base_C | base_A,
                'H' : base_T | base_G | base_A,
                'V' : base_T | base_G | base_C,
                'N' : base_T | base_G | base_A | base_C }


class SeqDeTagger:
    # Translation of bases for primer identification.
    def __init__(self, p5, p5s, p3, p3s, t5, t3, keep_primer=False):
        self.p5s = p5s
        self.p3s = p3s
        self.t5 = t5
        self.t3 = t3
        self.keep_primer = keep_primer

        # Translate primer sequence
        self.p5 = [trans_table[x] for x in p5.upper() if x in trans_table]
        self.p3 = [complement_trans_table[x] for x in p3.upper() if x in complement_trans_table]
        self.p3.reverse()
        self.p5len = float(len(self.p5))
        self.p3len = float(len(self.p3))


    # Function to detag sequence and return a dict()
    # with sequence and metadata

    def detag_seq(self, qualseq):
    
        result = DeTaggedSeq()
        seq_record = qualseq.get_seq()
        seq = seq_record.seq
        seq_str = str(seq)
        seq_list = [trans_table[x] if not x == 'N' else 0 for x in str(seq).upper() if x in trans_table]
        q=qualseq.get_qual()        

        p5_pos = -1
        accepted_t3 = set()
        if self.p5:
            p5_len = len(self.p5)
            for x in range(0, min(len(seq_str), 20000 + p5_len ) - p5_len):
                m = 0
                for i in range(p5_len):
                    if seq_list[x + i] & self.p5[i]:
                        m += 1
                if p5_len - m <= self.p5s:
                    p5_pos = x
                    break

            if p5_pos < 0:
                result.reversed = True
                seq = seq.reverse_complement()
                if q:
                    q.quals.reverse()
                seq_str = str(seq)
                seq_list = [trans_table[x] if not x == 'N' else 0 for x in str(seq).upper() if x in trans_table]

                for x in range(0, min(len(seq_str), 20000 + p5_len ) - p5_len):
                    m = 0
                    for i in range(p5_len):
                        if seq_list[x + i] & self.p5[i]:
                            m += 1
                    if p5_len - m <= self.p5s:
                        p5_pos = x
                        break
            if p5_pos < 0:
                raise ScataReadsError("no_primer5", "No 5' primer found")
            
        else:
            p5_pos=0
        
        if self.t5:
            tag_len = len(next(iter(self.t5)))
            tag_seq = seq_str[(p5_pos - tag_len):p5_pos]
            try:
                result.tag = self.t5[tag_seq]['name']
                accepted_t3 = self.t5[tag_seq]['mates']
            except KeyError:
                raise ScataReadsError("no_tag5", "No 5' tag found")
        else:
            result.tag = ""
        
        if not self.keep_primer:
            p5_pos += p5_len

        p3_pos = -1
        p3_matches = [ ]

        if self.p3:
            p3_len = len(self.p3)
            for x in range(len(seq) - p3_len - 1, p5_pos, -1):
                m=0
                for i in range(p3_len):
                    if seq_list[x + i] & self.p3[i]:
                        m += 1
                if p3_len - m <= self.p3s:
                    p3_pos = x
                    break
            if p3_pos < 0:
                raise ScataReadsError("no_primer3", "No 3' primer found")
    
            if self.t3:
                tag_len = len(next(iter(self.t3)))
                tag_seq = seq_str[(p3_pos + len(self.p3)):(p3_pos + len(self.p3) + tag_len)]
                tag_seq = str(Seq(tag_seq).reverse_complement())
                try:
                    if len(accepted_t3) and self.t3[tag_seq]['name'] not in accepted_t3:
                        raise ScataReadsError("chimeric_tag", "Chimeric tags")                    
                    result.tag += ("_" + self.t3[tag_seq]['name'])
                except KeyError:
                    raise ScataReadsError("no_tag3", "No 3' tag found")
            else:
                result.tag += ""

                if self.keep_primer:
                    p3_pos += p3_len

            seq_record.seq = seq[(p5_pos + len(self.p5)):p3_pos]
            result.seq_record = seq_record
            if q:
                q.quals = q.quals[(p5_pos + len(self.p5)):p3_pos]
                result.qual = q

        else:
            seq_record.seq = seq[(p5_pos + len(self.p5)):]
            result.seq_record = seq_record
            if q:
                q.quals = q.quals[(p5_pos + len(self.p5)):]
                result.qual = q

            
        return result

