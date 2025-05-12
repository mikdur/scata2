import gzip,pickle
from io import BytesIO
import time
from django.core.files import File
from scata2.models import ScataDataset, ScataErrorType
from scata2.backend.ReadHandler import Reads, ScataReadsError, ScataFileError
import scata2.backend
import django_q.tasks as q2



def check_dataset(pk):
    dataset = ScataDataset.objects.get(pk=pk)

    # Dictionary with all reads, with id as key
    #
    # {"read_id" : "read sequence"}
    #
    seqs = dict()

    # Dictionary of mappings between tag and readIDs
    # 
    #  {"tag_id":{"reads":[ .. ], "cnt": NN, "rev": NN}
    #
    tags = dict()

    try:
        file1 = gzip.open(dataset.file1.file.open(mode="rb"), mode="rt")
        if dataset.file2:
            file2 = gzip.open(dataset.file2.file.open(mode="rb"), mode="rt")
        else:
            file2 = None
    except gzip.BadGzipFile as e:
            dataset.validated = True
            dataset.is_valid = False
            dataset.progress = "Failed: not a gzipped file"
            dataset.save()
            return

    reads = Reads(file1=file1,
                  file2=file2,
                  file_type=dataset.file_types,
                  filtering=dataset.filter_method,
                  mean_min=dataset.mean_qual, min_qual=dataset.min_qual,
                  amplicon=dataset.amplicon,
                  kmer=dataset.kmer_size,
                  hsp=dataset.kmer_shared,
                  hsp_min=dataset.kmer_hsp_count)

    start_time = time.process_time()

    total_reads = 0
    good_reads = 0
    rev_reads = 0

    filter_results = dict()
    while True:
        try:
            total_reads += 1
            if total_reads % 10000 == 0:
                dataset.progress = ("Filtering, {t} reads done. {g} reads " +
                "accepted").format(g=good_reads, t=total_reads)
                dataset.save()
            read = next(reads)
            
            
            if read.tag in tags:
                tags[read.tag]["cnt"] += 1
                if read.reversed:
                    tags[read.tag]["rev"] += 1
                tags[read.tag]["reads"].append(read.seq_record.id)
            else:
                tags[read.tag] = {"cnt": 1,
                                  "rev": 1 if read.reversed else 0,
                                  "reads": [ read.seq_record.id ],
                                  }
        
            seqs[read.seq_record.id]=read.seq_record.seq.upper()
            good_reads += 1 
            if read.reversed:
                rev_reads +=1
            
        except ScataReadsError as e:
            if e.error in filter_results:
                filter_results[e.error]["cnt"] += 1
            else:
                filter_results[e.error] = {"msg":e.message,
                                           "cnt":1}
        except ScataFileError as e:
            dataset.validated = True
            dataset.is_valid = False
            dataset.progress = "Failed: " + e.message
            dataset.save()
            return
        except gzip.BadGzipFile as e:
            dataset.validated = True
            dataset.is_valid = False
            dataset.progress = "Failed: not a gzipped file/broken gzip file"
            dataset.save()
            return
        except StopIteration:
            break

    dataset.progress = "Finalising, {g}/{t} good reads".format(g=good_reads, t=total_reads)
    dataset.save()

    # Save data to files
    with BytesIO() as tag_file:
        with gzip.open(tag_file, mode="wb") as gz:
            pickle.dump(tags, gz)
        tag_file.seek(0)
        name = "tags_{id}".format(id=pk)
        dataset.tags.save(name, File(tag_file, name=name))

    with BytesIO() as seq_file:
        with gzip.open(seq_file, "wb") as gz:
            pickle.dump(seqs, gz)
        seq_file.seek(0)
        name = "seqs_{id}".format(id=pk)
        dataset.sequences.save(name, File(seq_file, name=name))
    


    dataset.seq_count = good_reads
    dataset.seq_total = total_reads
    dataset.seq_rev = rev_reads
    dataset.tag_count = len(tags)
    dataset.process_time = time.process_time() - start_time
    dataset.validated = True
    if good_reads > 0:
        dataset.is_valid = True
    dataset.progress = "Ready, {g}/{t} good reads".format(g=good_reads, t=total_reads)
    dataset.save()
    
    for e in filter_results.items():
        err = ScataErrorType()
        err.error = e[0]
        err.message = e[1]['msg']
        err.count = e[1]['cnt']
        err.dataset = dataset
        err.save()

    q2.async_task(scata2.backend.dataset_stats, pk,
                      task_name="dataset stats pk={id}".format(id=pk))    

