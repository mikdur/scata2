import gzip,pickle
from io import BytesIO, TextIOWrapper
import time
from django.core.files import File
from scata2.models import ScataReferenceSet, ScataRefsetErrorType
from scata2.backend.ReadHandler import Reads, ScataReadsError, ScataFileError

def check_refset(pk):
    refset = ScataReferenceSet.objects.get(pk=pk)

    # Dictionary with all reads, with id as key
    #
    # {"read_id" : "read sequence"}
    #
    seqs = dict()

    
    try:
        file1 = gzip.open(refset.refseq_file.file.open(mode="rb"), mode="rt")
        next(file1)
        file1 = gzip.open(refset.refseq_file.file.open(mode="rb"), mode="rt")

    except gzip.BadGzipFile as e:
            file1.close()
            file1 = None
            pass

    if not file1:
         file1 = TextIOWrapper(refset.refseq_file.file.open(mode="rb"))

    reads = Reads(file1=file1,
                  file_type="fasta",
                  filtering="fs",
                  amplicon=refset.amplicon,
                  keep_primer=False,
                  ignore_tags=True)
    

    start_time = time.process_time()

    total_reads = 0
    good_reads = 0

    filter_results = dict()

    while True:
        try:
            total_reads += 1
            if total_reads % 10000 == 0:
                refset.progress = ("Filtering, {t} sequences done. {g} sequences " +
                "accepted").format(g=good_reads, t=total_reads)
                refset.save()
            read = next(reads)
            
        
            seqs[read.seq_record.id]=read.seq_record.seq
            good_reads += 1     
            
        except ScataReadsError as e:
            if e.error in filter_results:
                filter_results[e.error]["cnt"] += 1
            else:
                filter_results[e.error] = {"msg":e.message,
                                           "cnt":1}
        except ScataFileError as e:
            refset.validated = True
            refset.is_valid = False
            refset.progress = "Failed: " + e.message
            refset.save()
            return
        except gzip.BadGzipFile as e:
            refset.validated = True
            refset.is_valid = False
            refset.progress = "Failed: not a gzipped file/broken gzip file"
            refset.save()
            return
        except StopIteration:
            break

    refset.progress = "Finalising, {g}/{t} good reads".format(g=good_reads, t=total_reads)
    refset.save()

    with BytesIO() as seq_file:
        with gzip.open(seq_file, "wb") as gz:
            pickle.dump(seqs, gz)
        seq_file.seek(0)
        name = "refs_{id}".format(id=pk)
        refset.sequences.save(name, File(seq_file, name=name))
    


    refset.seq_count = good_reads
    refset.seq_total = total_reads
    refset.process_time = time.process_time() - start_time
    refset.validated = True
    if good_reads > 0:
        refset.is_valid = True
    refset.progress = "Ready, {g}/{t} good reads".format(g=good_reads, t=total_reads)
    refset.save()

    print(filter_results)
    for e in filter_results.items():
        err = ScataRefsetErrorType()
        err.error = e[0]
        err.message = e[1]['msg']
        err.count = e[1]['cnt']
        err.refset = refset
        err.save()

