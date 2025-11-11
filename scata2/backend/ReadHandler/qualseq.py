
class QualSeq:
    def __init__(self, seq_record, qual=None):
        self.seq_record = seq_record
        self.qual = qual

    def __len__ (self):
        return len(self.seq_record.seq)

    def get_seq(self):
        return self.seq_record
    
    def get_qual(self):
        return self.qual



class Qual:
    def __init__ (self, name, qual):
        self.name = name
        self.quals = qual
    

    def __getitem__ (self, item):
        return int(self.quals[item])

    def __repr__ (self):
        return "Qual(name = " + self.name + ", quals = [" + (", ".join([str(x) for x in self.quals])) + "]"


class QualFile:
    def __init__ (self, qualfile):
        self.qualfile = qualfile
        try:
            self.line = next(qualfile)
        except StopIteration:
            raise Exception("Format error in qualfile")
        self.eof=0

    def __next__(self):
        if self.eof == 1:
            raise StopIteration

        quals = [ ]
        
        if not self.line[0] == ">":
            raise Scat ("Format error in qualfile")

        seq_name = self.line[1:].split()[0]
        
        while True:
            try:
                self.line = next(self.qualfile)
                
                if self.line[0] == ">":
                    return Qual(seq_name, quals)
        
                quals += [int(x) for x in self.line.split()]
            except StopIteration:
                self.eof = 1
                return(Qual(seq_name, quals))
