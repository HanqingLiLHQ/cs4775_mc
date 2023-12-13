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
    total_count=0.0
    motifs_fetched=[]
    cluster_01 = ['MA0119.1', 'MA1643.2', 'MA1527.2', 'MA0434.2', 'MA1712.2']
    motifs_cluster01 = []
    for ma_id in cluster_01:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster01.append(motif)
    motifs_fetched.append(motifs_cluster01)
    total_count += len(motifs_cluster01)

    cluster_02 = ['MA2323.1', 'MA2327.1', 'MA0113.4', 'MA0727.2', 'MA0007.4']
    motifs_cluster02 = []
    for ma_id in cluster_02:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster02.append(motif)
    motifs_fetched.append(motifs_cluster02)
    total_count += len(motifs_cluster02)

    cluster_03 = ['MA0813.1', 'MA0815.1', 'MA0872.1', 'MA0154.5', 'MA1637.2', 'MA0810.2','MA0811.2','MA0524.3','MA2122.1','MA1604.2']
    motifs_cluster03 = []
    for ma_id in cluster_03:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster03.append(motif)
    motifs_fetched.append(motifs_cluster03)
    total_count += len(motifs_cluster03)

    cluster_04 = ['MA0258.2', 'MA0112.4', 'MA0066.2']
    motifs_cluster04 = []
    for ma_id in cluster_04:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster04.append(motif)
    motifs_fetched.append(motifs_cluster04)
    total_count += len(motifs_cluster04)

    cluster_05 = ['MA0525.2', 'MA0106.3', 'MA0861.2']
    motifs_cluster05 = []
    for ma_id in cluster_05:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster05.append(motif)
    motifs_fetched.append(motifs_cluster05)
    total_count += len(motifs_cluster05)

    cluster_06 = ['MA0159.1', 'MA0730.1', 'MA0858.1', 'MA1149.2']
    motifs_cluster06 = []
    for ma_id in cluster_06:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster06.append(motif)
    motifs_fetched.append(motifs_cluster06)
    total_count += len(motifs_cluster06)

    



    # clusters=[
    # cluster_31, cluster_32, cluster_33, cluster_34, cluster_35, cluster_36, cluster_37, cluster_38, cluster_39, cluster_40,
    # cluster_41, cluster_42, cluster_43, cluster_44, cluster_45, cluster_46, cluster_47, cluster_48, cluster_49, cluster_50,
    # cluster_51, cluster_52, cluster_53, cluster_54, cluster_55, cluster_56, cluster_57, cluster_58, cluster_59, cluster_60,
    # cluster_61, cluster_62, cluster_63, cluster_64, cluster_65
    # ]
    # i=31
    # for cluster in clusters:
    #   cluster_tem=[]
    #   for ma_id in cluster:
    #     motif = jdb.fetch_motif_by_id(ma_id)
    #     cluster_tem.append(motif)
    #   total_count=total_count+len(cluster_tem)
    #   motifs_fetched.append(cluster_tem)
    # self.total_cluster_count = len(motifs_fetched)
    # print(len(motifs_fetched))
    # print(total_count)
    # self.total_motif_count=total_count
    # self.dbs=motifs_fetched
    
    self.total_cluster_count=6
    self.total_motif_count=total_count
    print(str(total_count)+"motifs have been fetched.")
    self.dbs=motifs_fetched
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
    correct_clustering  = []
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


#db = DNABindingMotifs()
#sth1 = db.mmm[list(db.mmm.keys())[0]]
#sth2 = db.mmm[list(db.mmm.keys())[0]]
#print(sth1.pssm)
# sth1.pseudocounts = motifs.jaspar.calculate_pseudocounts(sth1)
# sth2.pseudocounts = {'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
# print(sth1.pssm)
# print(sth1.pssm.dist_pearson(sth2.pssm))
