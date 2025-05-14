import gzip,pickle
from io import BytesIO
import time
from django.core.files import File
from scata2.models import ScataDataset, ScataTagStat

from sklearn.feature_extraction import FeatureHasher
from sklearn import preprocessing
from sklearn.decomposition import PCA


# Kmer-length for PCA plot
KMER=5
KMER_MIN_SEQS=200


def dataset_stats(pk):
    dataset = ScataDataset.objects.get(pk=pk)

    # Dictionary with all reads, with id as key
    #
    # {"read_id" : "read sequence"}
    #
    seqs = dict()

    with dataset.sequences.file.open(mode="rb") as f:
        with gzip.open(f, mode="rb") as gz:
            seqs = pickle.load(gz)



    # Dictionary of mappings between tag and readIDs
    # 
    #  {"tag_id":{"reads":[ .. ], "cnt": NN, "rev": NN}
    #
    tags = dict()

    with dataset.tags.file.open(mode="rb") as f:
        with gzip.open(f, mode="rb") as gz:
            tags = pickle.load(gz)

    # Make static list of tags (to ensure order is stable)

    tag_list = tags.keys()
    tag_objects = [] 

    kmer_list = []
    pca_objects = []

    # Calculate stats per tag
    for t in tag_list:
       
        tag = ScataTagStat()
        tag.dataset=dataset
        tag.count = tags[t]['cnt']
        tag.reversed = tags[t]['rev']

        # Scan all sequences in tag and calculate
        # min/max/mean length
        # min/max/mean gc
        # kmer composition
       
        lens = []
        gcs = []
        kmers = dict()

        for s in tags[t]['reads']:
            seq = str(seqs[s])
            lens.append(len(seq))
            gc = 0
            for i in range(len(seq)):
                if seq[i] in ([ "G", "C" ]):
                    gc += 1
                # Get kmers
                if i < len(seq) - KMER:
                    kmer = seq[i:i+KMER]
                    kmers[kmer] = kmers.get(kmer, 0) + 1
            gcs.append(gc / len(seq))


        tag.tag = t
        tag.min_len = min(lens)
        tag.max_len = max(lens)
        tag.mean_len = sum(lens) / len(lens)
        tag.min_gc = min(gcs)
        tag.max_gc = max(gcs)
        tag.mean_gc = sum(gcs) / len(gcs)
        tag_objects.append(tag)

        if tag.count > KMER_MIN_SEQS:
            kmer_list.append(kmers)
            pca_objects.append(tag)

    
    hasher = FeatureHasher(n_features=4**KMER, input_type="dict")

    f = hasher.transform(kmer_list)
    f = preprocessing.normalize(f, copy=False)
    pca = PCA(n_components=3)
    pca.fit(f)

    pca.explained_variance_ratio_
    eigen_vectors = pca.transform(f)

    dataset.pc1_exp=float(pca.explained_variance_ratio_[0])
    dataset.pc2_exp=float(pca.explained_variance_ratio_[1])
    dataset.pc3_exp=float(pca.explained_variance_ratio_[2])
    dataset.save()
    
    for i in range(len(pca_objects)):
        pca_objects[i].pc1=float(eigen_vectors[i,0])
        pca_objects[i].pc2=float(eigen_vectors[i,1])
        pca_objects[i].pc3=float(eigen_vectors[i,2])
        pca_objects[i].in_pca = True

    # Persist all data to database

    for i in range(len(tag_objects)):
        tag_objects[i].save()
    



    