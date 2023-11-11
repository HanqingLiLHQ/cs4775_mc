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
    cluster_1 = [
    "MA0268.1", "MA0306.1", "MA0337.1", "MA0338.1", "MA0339.1", "MA0341.1",
    "MA0342.1", "MA0364.1", "MA0366.1", "MA0372.1", "MA0423.1", "MA0431.1",
    "MA0436.1", "MA0441.1", "MA1431.1", "MA1432.1", "MA1433.1", "MA1434.1"
    ]
    motifs_cluster1 = []
    for ma_id in cluster_1:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster1.append(motif)
    motifs_fetched.append(motifs_cluster1)
    total_count=total_count+len(motifs_cluster1)

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

    cluster_6 = [
    'MA0275.1', 'MA0280.1', 'MA0282.1', 'MA0292.1', 'MA0294.1', 'MA0311.1', 'MA0325.1', 'MA0348.1', 'MA0354.1', 'MA0360.1', 'MA0367.1',
    'MA0380.1', 'MA0391.1', 'MA0392.1', 'MA0420.1', 'MA0422.1', 'MA0424.1', 'MA0428.1', 'MA0430.1', 'MA0432.1', 'MA0437.1', 'MA0438.1', 'MA0439.1'
    ]
    motifs_cluster6 = []
    for ma_id in cluster_6:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster6.append(motif)
    motifs_fetched.append(motifs_cluster6)
    total_count=total_count+len(motifs_cluster6)

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

    cluster_11 = ['MA0267', 'MA0402']
    motifs_cluster11 = []
    for ma_id in cluster_11:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster11.append(motif)
    motifs_fetched.append(motifs_cluster11)
    total_count=total_count+len(motifs_cluster11)

    cluster_12 = ['MA0289.1', 'MA0300.1', 'MA0307.1', 'MA0309.1', 'MA0312.2']
    motifs_cluster12 = []
    for ma_id in cluster_12:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster12.append(motif)
    motifs_fetched.append(motifs_cluster12)
    total_count=total_count+len(motifs_cluster12)

    cluster_13 = ['MA0355.1', 'MA0385.1']
    motifs_cluster13 = []
    for ma_id in cluster_13:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster13.append(motif)
    motifs_fetched.append(motifs_cluster13)
    total_count=total_count+len(motifs_cluster13)


    cluster_14 = ['MA0304.1', 'MA0305.1']
    motifs_cluster14 = []
    for ma_id in cluster_14:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster14.append(motif)
    motifs_fetched.append(motifs_cluster14)
    total_count=total_count+len(motifs_cluster14)


    cluster_15 = ['MA0344.1', 'MA0358.1', 'MA0421.1']
    motifs_cluster15 = []
    for ma_id in cluster_15:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster15.append(motif)
    motifs_fetched.append(motifs_cluster15)
    total_count=total_count+len(motifs_cluster15)

    cluster_16 = ['MA0284.2', 'MA0382.2', 'MA0415.1', 'MA0418.1']
    motifs_cluster16 = []
    for ma_id in cluster_16:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster16.append(motif)
    motifs_fetched.append(motifs_cluster16)
    total_count=total_count+len(motifs_cluster16)
    
    cluster_17 = ['MA0279.2', 'MA0297.1', 'MA0317.1', 'MA0419.1', 'MA0929.1']
    motifs_cluster17 = []
    for ma_id in cluster_17:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster17.append(motif)
    motifs_fetched.append(motifs_cluster17)
    total_count=total_count+len(motifs_cluster17)

    cluster_18 = [
    'MA0283.1', 'MA0320.1', 'MA0353.1', 'MA0374.1',
    'MA0375.1', 'MA0394.1', 'MA0399.1', 'MA0410.1', 'MA1435.1'
    ]
    motifs_cluster18 = []
    for ma_id in cluster_18:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster18.append(motif)
    motifs_fetched.append(motifs_cluster18)
    total_count=total_count+len(motifs_cluster18)

    cluster_19 = ['MA0308.1', 'MA0400.1']
    motifs_cluster19 = []
    for ma_id in cluster_19:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster19.append(motif)
    motifs_fetched.append(motifs_cluster19)
    total_count=total_count+len(motifs_cluster19)

    cluster_20 = ['MA0291.1', 'MA0405.1']
    motifs_cluster20 = []
    for ma_id in cluster_20:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster20.append(motif)
    motifs_fetched.append(motifs_cluster20)
    total_count=total_count+len(motifs_cluster20)

    cluster_21 =  ['MA0276.1', 'MA0324.1', 'MA0331.1']
    motifs_cluster21 = []
    for ma_id in cluster_21:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster21.append(motif)
    motifs_fetched.append(motifs_cluster21)
    total_count=total_count+len(motifs_cluster21)


    cluster_22 = ['MA0286.1', 'MA0416.1']
    motifs_cluster22 = []
    for ma_id in cluster_22:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster22.append(motif)
    motifs_fetched.append(motifs_cluster22)
    total_count=total_count+len(motifs_cluster22)

    cluster_23 = ['MA0336.1', 'MA0377.1']
    motifs_cluster23 = []
    for ma_id in cluster_23:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster23.append(motif)
    motifs_fetched.append(motifs_cluster23)
    total_count=total_count+len(motifs_cluster23)

    cluster_24 =  ['MA0313.1', 'MA0368.1']
    motifs_cluster24 = []
    for ma_id in cluster_24:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster24.append(motif)
    motifs_fetched.append(motifs_cluster24)
    total_count=total_count+len(motifs_cluster24)

    cluster_25 = ['MA0271.1', 'MA0272.1', 'MA0303.2']
    motifs_cluster25 = []
    for ma_id in cluster_25:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster25.append(motif)
    motifs_fetched.append(motifs_cluster25)
    total_count=total_count+len(motifs_cluster25)

    cluster_26 = ['MA0281.2', 'MA0376.1']
    motifs_cluster26 = []
    for ma_id in cluster_26:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster26.append(motif)
    motifs_fetched.append(motifs_cluster26)
    total_count=total_count+len(motifs_cluster26)

    cluster_27 = ['MA0274.1', 'MA0356.1', 'MA0387.1', 'MA0426.1', 'MA0433.1']
    motifs_cluster27 = []
    for ma_id in cluster_27:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster27.append(motif)
    motifs_fetched.append(motifs_cluster27)
    total_count=total_count+len(motifs_cluster27)

    cluster_28 = ['MA0362.1', 'MA0414.1']
    motifs_cluster28 = []
    for ma_id in cluster_28:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster28.append(motif)
    motifs_fetched.append(motifs_cluster28)
    total_count=total_count+len(motifs_cluster28)

    cluster_29 = ['MA0349.1', 'MA0404.1']
    motifs_cluster29 = []
    for ma_id in cluster_29:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster29.append(motif)
    motifs_fetched.append(motifs_cluster29)
    total_count=total_count+len(motifs_cluster29)

    cluster_30 = ['MA0319.1', 'MA0393.1']
    motifs_cluster30 = []
    for ma_id in cluster_30:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster30.append(motif)
    motifs_fetched.append(motifs_cluster30)
    total_count=total_count+len(motifs_cluster30)

    cluster_31 = ['MA0288.1', 'MA0327.1']
    motifs_cluster31 = []
    for ma_id in cluster_31:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster31.append(motif)
    motifs_fetched.append(motifs_cluster31)
    total_count=total_count+len(motifs_cluster31)

    cluster_32 = ['MA0277.1', 'MA0388.1']
    motifs_cluster32 = []
    for ma_id in cluster_32:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster32.append(motif)
    motifs_fetched.append(motifs_cluster32)
    total_count=total_count+len(motifs_cluster32)

    cluster_33 = ['MA0378.1', 'MA0390.1']
    motifs_cluster33 = []
    for ma_id in cluster_33:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster33.append(motif)
    motifs_fetched.append(motifs_cluster33)
    total_count=total_count+len(motifs_cluster33)

    cluster_34 = ['MA0314.2', 'MA0316.1']
    motifs_cluster34 = []
    for ma_id in cluster_34:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster34.append(motif)
    motifs_fetched.append(motifs_cluster34)
    total_count=total_count+len(motifs_cluster34)

    cluster_35 = ['MA0323.1', 'MA0425.1']
    motifs_cluster35 = []
    for ma_id in cluster_35:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster35.append(motif)
    motifs_fetched.append(motifs_cluster35)
    total_count=total_count+len(motifs_cluster35)

    cluster_36 = ['MA0398.1']
    motifs_cluster36 = []
    for ma_id in cluster_36:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster36.append(motif)
    motifs_fetched.append(motifs_cluster36)
    total_count=total_count+len(motifs_cluster36)

    cluster_37 = ['MA0347.2', 'MA0403.2']
    motifs_cluster37 = []
    for ma_id in cluster_37:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster37.append(motif)
    motifs_fetched.append(motifs_cluster37)
    total_count=total_count+len(motifs_cluster37)

    cluster_38 = ['MA0340.1', 'MA0417.1']
    motifs_cluster38 = []
    for ma_id in cluster_38:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster38.append(motif)
    motifs_fetched.append(motifs_cluster38)
    total_count=total_count+len(motifs_cluster38)

    cluster_39 = ['MA0278.1']
    motifs_cluster39 = []
    for ma_id in cluster_39:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster39.append(motif)
    motifs_fetched.append(motifs_cluster39)
    total_count=total_count+len(motifs_cluster39)

    cluster_40 = ['MA0363.2', 'MA0412.2']
    motifs_cluster40 = []
    for ma_id in cluster_40:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster40.append(motif)
    motifs_fetched.append(motifs_cluster40)
    total_count=total_count+len(motifs_cluster40)

    cluster_41 = ['MA0365.1']
    motifs_cluster41 = []
    for ma_id in cluster_41:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster41.append(motif)
    motifs_fetched.append(motifs_cluster41)
    total_count=total_count+len(motifs_cluster41)

    cluster_42 = ['MA0266.1']
    motifs_cluster42 = []
    for ma_id in cluster_42:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster42.append(motif)
    motifs_fetched.append(motifs_cluster42)
    total_count=total_count+len(motifs_cluster42)

    cluster_43 = ['MA0270.1']
    motifs_cluster43 = []
    for ma_id in cluster_43:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster43.append(motif)
    motifs_fetched.append(motifs_cluster43)
    total_count=total_count+len(motifs_cluster43)

    cluster_44 = ['MA0296.1']
    motifs_cluster44 = []
    for ma_id in cluster_44:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster44.append(motif)
    motifs_fetched.append(motifs_cluster44)
    total_count=total_count+len(motifs_cluster44)

    cluster_45 = ['MA0343.1']
    motifs_cluster45 = []
    for ma_id in cluster_45:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster45.append(motif)
    motifs_fetched.append(motifs_cluster45)
    total_count=total_count+len(motifs_cluster45)

    cluster_46 = ['MA0386.1']
    motifs_cluster46 = []
    for ma_id in cluster_46:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster46.append(motif)
    motifs_fetched.append(motifs_cluster46)
    total_count=total_count+len(motifs_cluster46)

    cluster_47 = ['MA0369.1']
    motifs_cluster47 = []
    for ma_id in cluster_47:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster47.append(motif)
    motifs_fetched.append(motifs_cluster47)
    total_count=total_count+len(motifs_cluster47)

    cluster_48 = ['MA0287.1']
    motifs_cluster48 = []
    for ma_id in cluster_48:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster48.append(motif)
    motifs_fetched.append(motifs_cluster48)
    total_count=total_count+len(motifs_cluster48)

    cluster_49 = ['MA0371.1']
    motifs_cluster49 = []
    for ma_id in cluster_49:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster49.append(motif)
    motifs_fetched.append(motifs_cluster49)
    total_count=total_count+len(motifs_cluster49)

    cluster_50 = ['MA0440.1']
    motifs_cluster50 = []
    for ma_id in cluster_50:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster50.append(motif)
    motifs_fetched.append(motifs_cluster50)
    total_count=total_count+len(motifs_cluster50)

    cluster_51 = ['MA0370.1']
    motifs_cluster51 = []
    for ma_id in cluster_51:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster51.append(motif)
    motifs_fetched.append(motifs_cluster51)
    total_count=total_count+len(motifs_cluster51)

    cluster_52 = ['MA0290.1']
    motifs_cluster52 = []
    for ma_id in cluster_52:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster52.append(motif)
    motifs_fetched.append(motifs_cluster52)
    total_count=total_count+len(motifs_cluster52)

    cluster_53 = ['MA0379.1', 'MA0411.1']
    motifs_cluster53 = []
    for ma_id in cluster_53:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster53.append(motif)
    motifs_fetched.append(motifs_cluster53)
    total_count=total_count+len(motifs_cluster53)

    cluster_54 = ['MA0406.1']
    motifs_cluster54 = []
    for ma_id in cluster_54:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster54.append(motif)
    motifs_fetched.append(motifs_cluster54)
    total_count=total_count+len(motifs_cluster54)

    cluster_55 = ['MA0273.1']
    motifs_cluster55 = []
    for ma_id in cluster_55:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster55.append(motif)
    motifs_fetched.append(motifs_cluster55)
    total_count=total_count+len(motifs_cluster55)

    cluster_56 = ['MA0352.2']
    motifs_cluster56 = []
    for ma_id in cluster_56:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster56.append(motif)
    motifs_fetched.append(motifs_cluster56)
    total_count=total_count+len(motifs_cluster56)

    cluster_57 = ['MA0326.1']
    motifs_cluster57 = []
    for ma_id in cluster_57:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster57.append(motif)
    motifs_fetched.append(motifs_cluster57)
    total_count=total_count+len(motifs_cluster57)

    cluster_58 = ['MA0413.1']
    motifs_cluster58 = []
    for ma_id in cluster_58:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster58.append(motif)
    motifs_fetched.append(motifs_cluster58)
    total_count=total_count+len(motifs_cluster58)

    cluster_59 = ['MA0435.1']
    motifs_cluster59 = []
    for ma_id in cluster_59:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster59.append(motif)
    motifs_fetched.append(motifs_cluster59)
    total_count=total_count+len(motifs_cluster59)

    cluster_60 = ['MA0359.2']
    motifs_cluster60 = []
    for ma_id in cluster_60:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster60.append(motif)
    motifs_fetched.append(motifs_cluster60)
    total_count=total_count+len(motifs_cluster60)

    cluster_61 = ['MA0269.1']
    motifs_cluster61 = []
    for ma_id in cluster_61:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster61.append(motif)
    motifs_fetched.append(motifs_cluster61)
    total_count=total_count+len(motifs_cluster61)

    cluster_62 = ['MA0395.1']
    motifs_cluster62 = []
    for ma_id in cluster_62:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster62.append(motif)
    motifs_fetched.append(motifs_cluster62)
    total_count=total_count+len(motifs_cluster62)

    cluster_63 = ['MA0265.2']
    motifs_cluster63 = []
    for ma_id in cluster_63:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster63.append(motif)
    motifs_fetched.append(motifs_cluster63)
    total_count=total_count+len(motifs_cluster63)

    cluster_64 = ['MA0299.1']
    motifs_cluster64 = []
    for ma_id in cluster_64:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster64.append(motif)
    motifs_fetched.append(motifs_cluster64)
    total_count=total_count+len(motifs_cluster64)

    cluster_65 = ['MA0407.1']
    motifs_cluster65 = []
    for ma_id in cluster_65:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster65.append(motif)
    motifs_fetched.append(motifs_cluster65)
    total_count=total_count+len(motifs_cluster65)

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


db = DNABindingMotifs()
sth1 = db.mmm[list(db.mmm.keys())[0]]
sth2 = db.mmm[list(db.mmm.keys())[0]]
print(sth1.pssm)
sth1.pseudocounts = motifs.jaspar.calculate_pseudocounts(sth1)
sth2.pseudocounts = {'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
print(sth1.pssm)
print(sth1.pssm.dist_pearson(sth2.pssm))
