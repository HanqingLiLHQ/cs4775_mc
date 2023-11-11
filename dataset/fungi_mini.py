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
    

    cluster_2 = ['MA0293.1', 'MA0301.1', 'MA0302.1', 'MA0389.1', 'MA0434.1', 'MA1437.1']
    motifs_cluster2 = []
    for ma_id in cluster_2:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster2.append(motif)
    motifs_fetched.append(motifs_cluster2)
    total_count=total_count+len(motifs_cluster2)
    
    cluster_3 = ['MA0384.1', 'MA0396.1', 'MA0397.1']
    motifs_cluster3 = []
    for ma_id in cluster_3:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster3.append(motif)
    motifs_fetched.append(motifs_cluster3)
    total_count=total_count+len(motifs_cluster3)

    cluster_4 = ['MA0318.1', 'MA0321.1', 'MA0322.1', 'MA0328.2', 'MA0408.1']
    motifs_cluster4 = []
    for ma_id in cluster_4:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster4.append(motif)
    motifs_fetched.append(motifs_cluster4)  
    total_count=total_count+len(motifs_cluster4)

    cluster_5 = ['MA0350.1', 'MA0351.1']
    motifs_cluster5 = []
    for ma_id in cluster_5:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster5.append(motif)
    motifs_fetched.append(motifs_cluster5)
    total_count=total_count+len(motifs_cluster5)

    

    cluster_7 = ['MA0310.1', 'MA0357.1', 'MA0409.1', 'MA1436.1']
    motifs_cluster7 = []
    for ma_id in cluster_7:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster7.append(motif)
    motifs_fetched.append(motifs_cluster7)
    total_count=total_count+len(motifs_cluster7)

    cluster_8 = ['MA0285.1', 'MA0361.1', 'MA0381.1', 'MA0429.1']
    motifs_cluster8 = []
    for ma_id in cluster_8:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster8.append(motif)
    motifs_fetched.append(motifs_cluster8)
    total_count=total_count+len(motifs_cluster8)

    cluster_9 = ['MA0295.1', 'MA0329.1', 'MA0330.1', 'MA0401.1']
    motifs_cluster9 = []
    for ma_id in cluster_9:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster9.append(motif)
    motifs_fetched.append(motifs_cluster9)
    total_count=total_count+len(motifs_cluster9)

    cluster_10 = ['MA0332.1', 'MA0333.1', 'MA0334.1', 'MA0373.1']
    motifs_cluster10 = []
    for ma_id in cluster_10:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster10.append(motif)
    motifs_fetched.append(motifs_cluster10)
    total_count=total_count+len(motifs_cluster10)

    
    

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
    
    self.total_cluster_count=65
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


# db = DNABindingMotifs()
# sth1 = db.mmm[list(db.mmm.keys())[0]]
# sth2 = db.mmm[list(db.mmm.keys())[0]]
# print(sth1.pssm)
# sth1.pseudocounts = motifs.jaspar.calculate_pseudocounts(sth1)
# sth2.pseudocounts = {'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
# print(sth1.pssm)
# print(sth1.pssm.dist_pearson(sth2.pssm))
