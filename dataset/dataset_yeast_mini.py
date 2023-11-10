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
        # fetch: yeast, Saccharomyces cerevisiae bZIP, bHLH
        # bHLH, bHSH, leucine zippers, homeodomain, and zinc finger factors are fetched.
        bHLH_s_cerevisiae_ABF1_1 = jdb.fetch_motif_by_id("MA0265.1")
        print("fetched first form of a bHLH motif ABF1")

        bHLH_s_cerevisiae_ABF1_2 = jdb.fetch_motif_by_id("MA0265.2")
        print("fetched second form of a bHLH motif ABF1")

        bHLH_db_s_cerevisiae = [
            bHLH_s_cerevisiae_ABF1_1, bHLH_s_cerevisiae_ABF1_2]

        CZF_s_cerevisiae_ACE2 = jdb.fetch_motif_by_id("MA0267.1")
        print("fetched a C2H2 zinc finger factors motif ACE2")

        CZF_s_cerevisiae_ADR1 = jdb.fetch_motif_by_id("MA0268.1")
        print("fetched a C2H2 zinc finger factors motif ADR1")

        CZF_s_cerevisiae = [CZF_s_cerevisiae_ACE2, CZF_s_cerevisiae_ADR1]

        bZIP_s_cerevisiae_ARR1 = jdb.fetch_motif_by_id("MA0274.1")
        print("fetched a bZIP motif ARR1")

        bZIP_s_cerevisiae_CAD1 = jdb.fetch_motif_by_id("MA0279.1")
        print("fetched a bZIP motif CAD1")

        bZIP_db_s_cerevisiae = [bZIP_s_cerevisiae_ARR1, bZIP_s_cerevisiae_CAD1]

        homeo_s_cerevisiae_HMRA1 = jdb.fetch_motif_by_id("MA0327.1")
        print("fetched a Homeo domain motif HMRA1")

        homeo_s_cerevisiae_HMRA2 = jdb.fetch_motif_by_id("MA0318.1")
        print("fetched a Homeo domain motif HMRA2")

        homeo_s_cerevisiae_MATALPHA2 = jdb.fetch_motif_by_id("MA0328.1")
        print("fetched a Homeo domain motif MATALPHA2")

        Homeo_db_s_cerevisiae = [homeo_s_cerevisiae_HMRA1,
                                 homeo_s_cerevisiae_HMRA2, homeo_s_cerevisiae_MATALPHA2]

        self.total_motif_count = len(bHLH_db_s_cerevisiae) + len(
            CZF_s_cerevisiae) + len(bZIP_db_s_cerevisiae) + len(Homeo_db_s_cerevisiae)
        self.total_cluster_count = 4
        print("fetched " + str(self.total_motif_count) + " motifs in total. ")
        self.dbs = [bHLH_db_s_cerevisiae, CZF_s_cerevisiae,
                    bZIP_db_s_cerevisiae, Homeo_db_s_cerevisiae]
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
# sth2.pseudocounts = {'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
# print(sth1.pssm)
# print(sth1.pssm.dist_pearson(sth2.pssm))
