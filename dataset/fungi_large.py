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
    cluster_01 = ['MA0310.2', 'MA0357.1', 'MA0376.2', 'MA0409.1', 'MA1436.2', 'MA0281.3', 'MA0321.1', 'MA0322.1', 'UN0089.2', 'UN0083.2', 'UN0294.2', 'UN0278.2', 'UN0293.1']
    motifs_cluster_01 = []
    for ma_id in cluster_01:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_01.append(motif)
    motifs_fetched.append(motifs_cluster_01)
    total_count=total_count+len(motifs_cluster_01)

    cluster_02 = ['MA0293.2', 'MA0302.2', 'MA0301.2', 'MA0389.2', 'MA0434.2', 'MA1437.2', 'MA0349.2', 'MA0404.1', 'UN0097.2']
    motifs_cluster_02 = []
    for ma_id in cluster_02:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_02.append(motif)
    motifs_fetched.append(motifs_cluster_02)
    total_count=total_count+len(motifs_cluster_02)

    cluster_03 = ['MA0381.2', 'MA0361.2', 'MA0429.1', 'MA0332.1', 'MA0334.2', 'MA0394.1', 'MA0395.2', 'MA0333.2', 'MA0285.2', 'MA0373.1', 'UN0088.1', 'UN0082.2', 'UN0081.2', 'UN0274.2', 'UN0094.2', 'UN0095.2', 'UN0103.2', 'UN0099.1']
    motifs_cluster_03 = []
    for ma_id in cluster_03:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_03.append(motif)
    motifs_fetched.append(motifs_cluster_03)
    total_count=total_count+len(motifs_cluster_03)

    cluster_04 = ['MA0296.2', 'MA0393.1', 'MA0317.2', 'MA0297.1', 'MA0929.2', 'UN0435.2']
    motifs_cluster_04 = []
    for ma_id in cluster_04:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_04.append(motif)
    motifs_fetched.append(motifs_cluster_04)
    total_count=total_count+len(motifs_cluster_04)

    cluster_05 = ['MA0379.1', 'MA0312.3', 'MA0307.1', 'MA0300.2', 'MA0289.1', 'MA0309.2']
    motifs_cluster_05 = []
    for ma_id in cluster_05:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_05.append(motif)
    motifs_fetched.append(motifs_cluster_05)
    total_count=total_count+len(motifs_cluster_05)

    cluster_06 = ['MA0278.2', 'MA0272.2', 'MA0303.3', 'MA0271.1', 'MA0295.2', 'UN0281.2']
    motifs_cluster_06 = []
    for ma_id in cluster_06:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_06.append(motif)
    motifs_fetched.append(motifs_cluster_06)
    total_count=total_count+len(motifs_cluster_06)

    cluster_07 = ['MA0414.2', 'MA0411.2', 'UN0096.1', 'UN0302.1']
    motifs_cluster_07 = []
    for ma_id in cluster_07:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_07.append(motif)
    motifs_fetched.append(motifs_cluster_07)
    total_count=total_count+len(motifs_cluster_07)

    cluster_08 = ['MA0291.1', 'MA0405.1', 'MA0308.2', 'MA0400.2', 'MA0437.1', 'MA0294.2', 'MA0367.2', 'MA0439.2', 'MA0354.2', 'MA0438.2', 'MA0428.2', 'MA0432.1', 'MA0392.1', 'MA0348.2', 'MA0422.2', 'MA0430.2', 'MA0325.2', 'UN0300.1', 'UN0283.2', 'UN0291.2', 'UN0295.1']
    motifs_cluster_08 = []
    for ma_id in cluster_08:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_08.append(motif)
    motifs_fetched.append(motifs_cluster_08)
    total_count=total_count+len(motifs_cluster_08)

    cluster_09 = ['MA0390.2', 'MA0378.2', 'MA0398.2']
    motifs_cluster_09 = []
    for ma_id in cluster_09:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_09.append(motif)
    motifs_fetched.append(motifs_cluster_09)
    total_count=total_count+len(motifs_cluster_09)

    cluster_10 = ['MA0304.1', 'MA0305.2', 'UN0091.2']
    motifs_cluster_10 = []
    for ma_id in cluster_10:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_10.append(motif)
    motifs_fetched.append(motifs_cluster_10)
    total_count=total_count+len(motifs_cluster_10)

    cluster_11 = ['MA0417.1', 'MA0340.1', 'UN0093.2']
    motifs_cluster_11 = []
    for ma_id in cluster_11:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_11.append(motif)
    motifs_fetched.append(motifs_cluster_11)
    total_count=total_count+len(motifs_cluster_11)

    cluster_12 = ['MA0362.2', 'MA0273.2', 'UN0286.1', 'UN0282.1', 'UN0287.1', 'UN0086.2', 'UN0290.1', 'UN0289.1', 'UN0292.2']
    motifs_cluster_12 = []
    for ma_id in cluster_12:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_12.append(motif)
    motifs_fetched.append(motifs_cluster_12)
    total_count=total_count+len(motifs_cluster_12)

    cluster_13 = ['MA0314.3', 'MA0313.1', 'MA0368.1']
    motifs_cluster_13 = []
    for ma_id in cluster_13:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_13.append(motif)
    motifs_fetched.append(motifs_cluster_13)
    total_count=total_count+len(motifs_cluster_13)

    cluster_14 = ['MA0269.2', 'MA0270.2']
    motifs_cluster_14 = []
    for ma_id in cluster_14:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_14.append(motif)
    motifs_fetched.append(motifs_cluster_14)
    total_count=total_count+len(motifs_cluster_14)

    cluster_15 = ['MA0352.3', 'MA0320.1', 'MA0353.1', 'MA1435.2', 'MA0344.1', 'MA0358.2', 'UN0301.2']
    motifs_cluster_15 = []
    for ma_id in cluster_15:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_15.append(motif)
    motifs_fetched.append(motifs_cluster_15)
    total_count=total_count+len(motifs_cluster_15)

    cluster_16 = ['MA0408.2', 'UN0297.2']
    motifs_cluster_16 = []
    for ma_id in cluster_16:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_16.append(motif)
    motifs_fetched.append(motifs_cluster_16)
    total_count=total_count+len(motifs_cluster_16)

    cluster_17 = ['MA0282.2', 'MA0391.1', 'MA0360.2', 'MA0311.1', 'MA0275.1', 'MA0420.2', 'MA0280.2', 'MA0380.2', 'MA0424.2', 'MA0292.2', 'UN0098.2', 'UN0102.1', 'UN0101.2', 'UN0298.2', 'UN0085.2', 'UN0280.2', 'UN0284.2', 'UN0277.1', 'UN0299.1', 'UN0285.2']
    motifs_cluster_17 = []
    for ma_id in cluster_17:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_17.append(motif)
    motifs_fetched.append(motifs_cluster_17)
    total_count=total_count+len(motifs_cluster_17)

    cluster_18 = ['MA0426.1', 'MA0433.2', 'MA0356.1', 'UN0296.2']
    motifs_cluster_18 = []
    for ma_id in cluster_18:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_18.append(motif)
    motifs_fetched.append(motifs_cluster_18)
    total_count=total_count+len(motifs_cluster_18)

    cluster_19 = ['MA0418.2', 'MA0382.3', 'MA0284.3', 'MA0415.2', 'MA0286.2', 'MA0416.2']
    motifs_cluster_19 = []
    for ma_id in cluster_19:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_19.append(motif)
    motifs_fetched.append(motifs_cluster_19)
    total_count=total_count+len(motifs_cluster_19)

    cluster_20 = ['MA0318.1', 'MA0328.2', 'MA0363.3', 'MA0421.1']
    motifs_cluster_20 = []
    for ma_id in cluster_20:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_20.append(motif)
    motifs_fetched.append(motifs_cluster_20)
    total_count=total_count+len(motifs_cluster_20)

    cluster_21 = ['MA0401.1', 'MA0375.2', 'MA0330.1', 'MA0329.2', 'MA0374.2', 'UN0084.1']
    motifs_cluster_21 = []
    for ma_id in cluster_21:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_21.append(motif)
    motifs_fetched.append(motifs_cluster_21)
    total_count=total_count+len(motifs_cluster_21)

    cluster_22 = ['MA0384.2', 'MA0396.2', 'MA0397.2']
    motifs_cluster_22 = []
    for ma_id in cluster_22:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_22.append(motif)
    motifs_fetched.append(motifs_cluster_22)
    total_count=total_count+len(motifs_cluster_22)

    cluster_23 = ['MA0343.2', 'MA0288.1', 'MA0327.1']
    motifs_cluster_23 = []
    for ma_id in cluster_23:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_23.append(motif)
    motifs_fetched.append(motifs_cluster_23)
    total_count=total_count+len(motifs_cluster_23)

    cluster_24 = ['MA0319.2', 'MA0336.2', 'MA0377.2']
    motifs_cluster_24 = []
    for ma_id in cluster_24:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_24.append(motif)
    motifs_fetched.append(motifs_cluster_24)
    total_count=total_count+len(motifs_cluster_24)

    cluster_25 = ['MA0306.2', 'MA0364.1', 'MA0423.2', 'MA0372.2', 'MA1433.2', 'MA0431.2', 'MA0441.2', 'MA0436.2', 'MA0268.2', 'MA0337.2', 'MA0338.2', 'MA0339.2', 'MA0425.2', 'MA1431.2', 'MA1432.2', 'MA0366.1', 'MA0341.1', 'MA0342.1', 'UN0100.2']
    motifs_cluster_25 = []
    for ma_id in cluster_25:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_25.append(motif)
    motifs_fetched.append(motifs_cluster_25)
    total_count=total_count+len(motifs_cluster_25)

    cluster_26 = ['MA0279.3', 'MA0419.1', 'UN0090.2', 'UN0288.1']
    motifs_cluster_26 = []
    for ma_id in cluster_26:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_26.append(motif)
    motifs_fetched.append(motifs_cluster_26)
    total_count=total_count+len(motifs_cluster_26)

    cluster_27 = ['MA0283.1', 'MA0399.1', 'MA0410.2', 'UN0279.1']
    motifs_cluster_27 = []
    for ma_id in cluster_27:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_27.append(motif)
    motifs_fetched.append(motifs_cluster_27)
    total_count=total_count+len(motifs_cluster_27)

    cluster_28 = ['MA0355.2', 'MA0385.2']
    motifs_cluster_28 = []
    for ma_id in cluster_28:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_28.append(motif)
    motifs_fetched.append(motifs_cluster_28)
    total_count=total_count+len(motifs_cluster_28)

    cluster_29 = ['MA0267.2', 'MA0402.2']
    motifs_cluster_29 = []
    for ma_id in cluster_29:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_29.append(motif)
    motifs_fetched.append(motifs_cluster_29)
    total_count=total_count+len(motifs_cluster_29)

    cluster_30 = ['MA0276.1', 'MA0331.1']
    motifs_cluster_30 = []
    for ma_id in cluster_30:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_30.append(motif)
    motifs_fetched.append(motifs_cluster_30)
    total_count=total_count+len(motifs_cluster_30)

    cluster_31 = ['MA0350.2', 'MA0351.2']
    motifs_cluster_31 = []
    for ma_id in cluster_31:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_31.append(motif)
    motifs_fetched.append(motifs_cluster_31)
    total_count=total_count+len(motifs_cluster_31)

    cluster_32 = ['MA0265.3']
    motifs_cluster_32 = []
    for ma_id in cluster_32:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_32.append(motif)
    motifs_fetched.append(motifs_cluster_32)
    total_count=total_count+len(motifs_cluster_32)

    cluster_33 = ['MA0266.2']
    motifs_cluster_33 = []
    for ma_id in cluster_33:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_33.append(motif)
    motifs_fetched.append(motifs_cluster_33)
    total_count=total_count+len(motifs_cluster_33)

    cluster_34 = ['MA0274.1']
    motifs_cluster_34 = []
    for ma_id in cluster_34:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_34.append(motif)
    motifs_fetched.append(motifs_cluster_34)
    total_count=total_count+len(motifs_cluster_34)

    cluster_35 = ['MA0277.1']
    motifs_cluster_35 = []
    for ma_id in cluster_35:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_35.append(motif)
    motifs_fetched.append(motifs_cluster_35)
    total_count=total_count+len(motifs_cluster_35)

    cluster_36 = ['MA0287.1']
    motifs_cluster_36 = []
    for ma_id in cluster_36:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_36.append(motif)
    motifs_fetched.append(motifs_cluster_36)
    total_count=total_count+len(motifs_cluster_36)

    cluster_37 = ['MA0290.1']
    motifs_cluster_37 = []
    for ma_id in cluster_37:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_37.append(motif)
    motifs_fetched.append(motifs_cluster_37)
    total_count=total_count+len(motifs_cluster_37)

    cluster_38 = ['MA0299.1']
    motifs_cluster_38 = []
    for ma_id in cluster_38:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_38.append(motif)
    motifs_fetched.append(motifs_cluster_38)
    total_count=total_count+len(motifs_cluster_38)

    cluster_39 = ['MA0316.2']
    motifs_cluster_39 = []
    for ma_id in cluster_39:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_39.append(motif)
    motifs_fetched.append(motifs_cluster_39)
    total_count=total_count+len(motifs_cluster_39)

    cluster_40 = ['MA0323.1']
    motifs_cluster_40 = []
    for ma_id in cluster_40:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_40.append(motif)
    motifs_fetched.append(motifs_cluster_40)
    total_count=total_count+len(motifs_cluster_40)

    cluster_41 = ['MA0324.1']
    motifs_cluster_41 = []
    for ma_id in cluster_41:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_41.append(motif)
    motifs_fetched.append(motifs_cluster_41)
    total_count=total_count+len(motifs_cluster_41)

    cluster_42 = ['MA0326.1']
    motifs_cluster_42 = []
    for ma_id in cluster_42:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_42.append(motif)
    motifs_fetched.append(motifs_cluster_42)
    total_count=total_count+len(motifs_cluster_42)

    cluster_43 = ['MA0347.3']
    motifs_cluster_43 = []
    for ma_id in cluster_43:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_43.append(motif)
    motifs_fetched.append(motifs_cluster_43)
    total_count=total_count+len(motifs_cluster_43)

    cluster_44 = ['MA0359.3']
    motifs_cluster_44 = []
    for ma_id in cluster_44:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_44.append(motif)
    motifs_fetched.append(motifs_cluster_44)
    total_count=total_count+len(motifs_cluster_44)

    cluster_45 = ['MA0365.2']
    motifs_cluster_45 = []
    for ma_id in cluster_45:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_45.append(motif)
    motifs_fetched.append(motifs_cluster_45)
    total_count=total_count+len(motifs_cluster_45)

    cluster_46 = ['MA0369.2']
    motifs_cluster_46 = []
    for ma_id in cluster_46:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_46.append(motif)
    motifs_fetched.append(motifs_cluster_46)
    total_count=total_count+len(motifs_cluster_46)

    cluster_47 = ['MA0370.1']
    motifs_cluster_47 = []
    for ma_id in cluster_47:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_47.append(motif)
    motifs_fetched.append(motifs_cluster_47)
    total_count=total_count+len(motifs_cluster_47)

    cluster_48 = ['MA0371.1']
    motifs_cluster_48 = []
    for ma_id in cluster_48:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_48.append(motif)
    motifs_fetched.append(motifs_cluster_48)
    total_count=total_count+len(motifs_cluster_48)

    cluster_49 = ['MA0386.2']
    motifs_cluster_49 = []
    for ma_id in cluster_49:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_49.append(motif)
    motifs_fetched.append(motifs_cluster_49)
    total_count=total_count+len(motifs_cluster_49)

    cluster_50 = ['MA0387.1']
    motifs_cluster_50 = []
    for ma_id in cluster_50:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_50.append(motif)
    motifs_fetched.append(motifs_cluster_50)
    total_count=total_count+len(motifs_cluster_50)

    cluster_51 = ['MA0388.1']
    motifs_cluster_51 = []
    for ma_id in cluster_51:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_51.append(motif)
    motifs_fetched.append(motifs_cluster_51)
    total_count=total_count+len(motifs_cluster_51)

    cluster_52 = ['MA0403.3']
    motifs_cluster_52 = []
    for ma_id in cluster_52:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_52.append(motif)
    motifs_fetched.append(motifs_cluster_52)
    total_count=total_count+len(motifs_cluster_52)

    cluster_53 = ['MA0406.2']
    motifs_cluster_53 = []
    for ma_id in cluster_53:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_53.append(motif)
    motifs_fetched.append(motifs_cluster_53)
    total_count=total_count+len(motifs_cluster_53)

    cluster_54 = ['MA0407.1']
    motifs_cluster_54 = []
    for ma_id in cluster_54:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_54.append(motif)
    motifs_fetched.append(motifs_cluster_54)
    total_count=total_count+len(motifs_cluster_54)

    cluster_55 = ['MA0412.3']
    motifs_cluster_55 = []
    for ma_id in cluster_55:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_55.append(motif)
    motifs_fetched.append(motifs_cluster_55)
    total_count=total_count+len(motifs_cluster_55)

    cluster_56 = ['MA0413.2']
    motifs_cluster_56 = []
    for ma_id in cluster_56:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_56.append(motif)
    motifs_fetched.append(motifs_cluster_56)
    total_count=total_count+len(motifs_cluster_56)

    cluster_57 = ['MA0435.2']
    motifs_cluster_57 = []
    for ma_id in cluster_57:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_57.append(motif)
    motifs_fetched.append(motifs_cluster_57)
    total_count=total_count+len(motifs_cluster_57)

    cluster_58 = ['MA0440.1']
    motifs_cluster_58 = []
    for ma_id in cluster_58:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_58.append(motif)
    motifs_fetched.append(motifs_cluster_58)
    total_count=total_count+len(motifs_cluster_58)

    cluster_59 = ['MA1434.1']
    motifs_cluster_59 = []
    for ma_id in cluster_59:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_59.append(motif)
    motifs_fetched.append(motifs_cluster_59)
    total_count=total_count+len(motifs_cluster_59)

    cluster_60 = ['UN0087.2']
    motifs_cluster_60 = []
    for ma_id in cluster_60:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_60.append(motif)
    motifs_fetched.append(motifs_cluster_60)
    total_count=total_count+len(motifs_cluster_60)

    cluster_61 = ['UN0092.1']
    motifs_cluster_61 = []
    for ma_id in cluster_61:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_61.append(motif)
    motifs_fetched.append(motifs_cluster_61)
    total_count=total_count+len(motifs_cluster_61)



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
    
    self.total_cluster_count=61
    self.total_motif_count=total_count
    print(str(total_count)+"motifs have been fetched.")
    #print(motifs_fetched)
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
