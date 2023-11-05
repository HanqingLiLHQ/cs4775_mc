# generate the jaspar dataset for DNA binding motifs
# Assume that pyjaspar, biopython are already installed
from pyjaspar import jaspardb
from Bio import motifs
from Bio.motifs import Motif
import random


class DNABindingMotifs(object):
    """
    Dataset for DNABindingMotifs
    Attributes:
      dbs: list of jaspar databases involved in data
      total_motif_count: the amount of motifs we have for this specific clustering task
      total_cluster_count: the amount of clusters an algorithm should give 
      mmm: motif -> matrix, (motif matrix map), this would be used as input for most clustering algorithms
      cc: correct clustering results, denoted by a list containing several maps. Each map is a motif -> matrix map
        denoting an original cluster in jaspar dataset. This will be used for the evaluation rubrics.
    """

    def __init__(self):
        self._fetch_data()
        self._cluster_overlap_check()
        self._generate_mmm()

    def _fetch_data(self):
        # helper function to initialize the dataset
        # initialize jaspar database object
        jdb = jaspardb()
        # print the version of jaspar database you are using
        print("This program runs on database: " + jdb.release)
        # fetch: vertebrates, Homo sapiens with different DNA binding motifs
        # bHLH, bHSH, leucine zippers, homeodomain, and zinc finger factors are fetched.
        bHLH_db = jdb.fetch_motifs(collection="core", tax_group=['vertebrates'],
                                   species=9606, tf_class="Basic helix-loop-helix factors (bHLH)")
        print("fetched " + str(len(bHLH_db)) +
              " motifs from Basic helix-loop-helix factors (bHLH)")

        bHSH_db = jdb.fetch_motifs(collection="core", tax_group=['vertebrates'],
                                   species=9606, tf_class="Basic helix-span-helix factors (bHSH)")
        print("fetched " + str(len(bHSH_db)) +
              " motifs from Basic helix-span-helix factors (bHSH)")

        bZIP_db = jdb.fetch_motifs(collection="core", tax_group=['vertebrates'],
                                   species=9606, tf_class="Basic leucine zipper factors (bZIP)")
        print("fetched " + str(len(bZIP_db)) +
              " motifs from Basic leucine zipper factors (bZIP)")

        Homeo_db = jdb.fetch_motifs(collection="core", tax_group=['vertebrates'],
                                    species=9606, tf_class="Homeo domain factors")
        print("fetched " + str(len(Homeo_db)) +
              " motifs from Homeo domain factors")

        Zinc_db = jdb.fetch_motifs(collection="core", tax_group=['vertebrates'],
                                   species=9606, tf_class="C2H2 zinc finger factors")
        print("fetched " + str(len(Zinc_db)) +
              " motifs from C2H2 zinc finger factors")

        self.total_motif_count = len(
            bHLH_db) + len(bHSH_db) + len(bZIP_db) + len(Homeo_db) + len(Zinc_db)
        self.total_cluster_count = 5
        print("fetched " + str(self.total_motif_count) + " motifs in total. ")
        self.dbs = [bHLH_db, bHSH_db, bZIP_db, Homeo_db, Zinc_db]
        return

    def _cluster_overlap_check(self):
        # check whether there are overlapping motifs between different classes.
        ids = []
        for db in self.dbs:
            for motif in db:
                if motif.matrix_id not in ids:
                    ids.append(motif.matrix_id)
                else:
                    print("found overlap between clusters")
        if len(ids) == self.total_motif_count:
            print("there is no overlap between clusters. ")
        return

    """
  Generate the motif-matrix-map. Motif refers to the matrix_id on jaspar website. 
  matrix is the position frequency matrix. This would generate for all the motifs 
  in all databases, therefore could be used as an input to clustering algorithms. 
  It will also generate the "correct" clustering results according to jaspar dataset. 
  """

    def _generate_mmm(self):
        sequence = []
        correct_clustering = []
        for db in self.dbs:
            cluster = {}
            for motif in db:
                sequence.append((motif.matrix_id, Motif(counts=motif.counts)))
                cluster[motif.matrix_id] = Motif(counts=motif.counts)
            correct_clustering.append(cluster)
        # shuffle the order of the map by first creating a list and then shuffle it
        random.shuffle(sequence)
        mmm = {}
        for (id, m) in sequence:
            mmm[id] = m

        assert len(mmm.keys()) == self.total_motif_count
        assert len(correct_clustering) == len(self.dbs)
        self.mmm = mmm
        self.cc = correct_clustering
        return


# db = DNABindingMotifs()
# sth1 = db.mmm[list(db.mmm.keys())[0]]
# sth2 = db.mmm[list(db.mmm.keys())[0]]
# print(sth1.pssm)
# sth1.pseudocounts = motifs.jaspar.calculate_pseudocounts(sth1)
# print(sth1.pseudocounts, sth1.pwm)
# sth2.pseudocounts = {'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
# print(sth1.pssm)
# print(sth1.pssm.dist_pearson(sth2.pssm))
