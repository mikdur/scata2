# This file is a dummy clustering method to illustrate and test basic
# API for clustering methods.

def do_cluster(instance):
    print("Clustering {}".format(instance))
    instance.job.status = "Preparing"
    instance.job.save()
    instance.mean_len = 0
    instance.num_seqs = 0

    seq_iter = instance.get_seq_iterator()

    clusters = dict()

    instance.job.status = "Deduplicating"
    instance.job.save()

    for seq in seq_iter:
        instance.mean_len += len(seq[1])
        instance.num_seqs += 1
        clusters[seq[1]] = clusters.get(seq[1], 0) + 1
    instance.mean_len =  instance.mean_len / instance.num_seqs
    instance.num_clusters = len(clusters)
    instance.save()
    instance.job.status = "Done"
    instance.job.save()

    print(seq_iter.errors)
    print(instance.mean_len)
    print(instance.num_seqs)
    print(instance.num_clusters)


