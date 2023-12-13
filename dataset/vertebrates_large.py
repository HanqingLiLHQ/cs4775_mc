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
    cluster_001 = ['MA1592.2', 'MA1928.2', 'MA0835.3', 'MA0462.3', 'MA1634.2', 'MA1988.2', 'MA1101.3', 'MA0476.2', 'MA0489.3', 'MA1132.2', 'MA1142.2', 'MA1134.2', 'MA1144.2', 'MA1135.2', 'MA1138.2', 'MA1128.2', 'MA1141.2', 'MA1137.2', 'MA0099.4', 'MA1130.2', 'MA1633.2', 'MA0655.1', 'MA0477.3', 'MA0490.3', 'MA0491.3', 'MA0478.2', 'MA0841.2', 'MA1520.2', 'MA1521.2', 'MA0496.4', 'MA0089.3', 'MA0150.3', 'MA0501.2', 'MA0591.2', 'MA0659.4', 'MA0738.2', 'MA0739.2', 'MA0161.3', 'MA0670.2', 'MA0671.2', 'MA0498.3', 'MA0774.1', 'MA0775.2', 'MA0782.3', 'MA1114.2', 'UN0799.1', 'UN0673.2', 'UN0807.1', 'UN0570.2']
    motifs_cluster_001 = []
    for ma_id in cluster_001:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_001.append(motif)
    motifs_fetched.append(motifs_cluster_001)
    total_count=total_count+len(motifs_cluster_001)

    cluster_002 = ['MA1643.2', 'MA0119.1', 'MA1528.2', 'MA1527.2', 'UN0127.2', 'UN0666.2', 'UN0621.2', 'UN0600.2', 'UN0629.2']
    motifs_cluster_002 = []
    for ma_id in cluster_002:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_002.append(motif)
    motifs_fetched.append(motifs_cluster_002)
    total_count=total_count+len(motifs_cluster_002)

    cluster_003 = ['MA0041.3', 'MA0847.4', 'MA1487.3', 'MA1607.2', 'MA0846.2', 'MA0032.2', 'MA0845.1', 'MA0635.2', 'MA0877.4', 'MA0040.2', 'MA0047.4', 'MA2118.1', 'MA0850.1', 'MA0848.1', 'MA0849.1', 'MA0031.2', 'MA0480.3', 'MA0157.4', 'MA0852.3', 'MA0042.2', 'MA2117.1', 'MA1103.3', 'MA0481.4', 'MA1683.2', 'MA0033.2', 'MA1606.2', 'MA0851.2', 'MA0030.2', 'MA0593.2', 'MA0148.5', 'MA1489.1', 'MA0613.1', 'MA0614.1', 'MA0651.3', 'MA0873.2', 'MA0906.2', 'MA1506.2', 'MA0907.2', 'MA0905.2', 'MA0485.3', 'MA1503.2', 'MA0908.2', 'MA0911.2', 'MA0909.4', 'MA0465.3', 'MA0650.4', 'MA0901.3', 'MA0878.3', 'MA1473.2', 'MA0899.2', 'MA0913.3', 'UN0802.1', 'UN0612.2', 'UN0648.2']
    motifs_cluster_003 = []
    for ma_id in cluster_003:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_003.append(motif)
    motifs_fetched.append(motifs_cluster_003)
    total_count=total_count+len(motifs_cluster_003)

    cluster_004 = ['MA0646.2', 'MA0767.2', 'UN0211.2', 'UN0554.2', 'UN0548.2', 'UN0551.2', 'UN0547.2', 'UN0553.2']
    motifs_cluster_004 = []
    for ma_id in cluster_004:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_004.append(motif)
    motifs_fetched.append(motifs_cluster_004)
    total_count=total_count+len(motifs_cluster_004)

    cluster_005 = ['MA0836.3', 'MA0102.5', 'MA0838.1', 'MA0466.4', 'MA0837.3', 'MA0833.3', 'MA1636.2', 'MA0639.2', 'MA0843.2', 'MA0025.3', 'MA0043.4', 'MA1933.2', 'MA1939.2', 'UN0522.2', 'UN0504.2', 'UN0510.2', 'UN0505.2', 'UN0527.2', 'UN0523.2', 'UN0526.1']
    motifs_cluster_005 = []
    for ma_id in cluster_005:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_005.append(motif)
    motifs_fetched.append(motifs_cluster_005)
    total_count=total_count+len(motifs_cluster_005)

    cluster_006 = ['MA0822.1', 'MA1485.1', 'MA0783.1', 'MA1571.1', 'MA1572.1', 'MA0796.1', 'MA0797.1', 'MA0103.4', 'MA0820.2', 'MA1558.2', 'MA0522.4', 'MA1648.2', 'MA0499.3', 'MA1631.2', 'MA1559.2', 'MA0830.3', 'MA0745.3', 'MA1620.2', 'MA0100.4', 'MA0521.3', 'MA1635.2', 'MA1997.2', 'MA0633.3', 'MA1993.2', 'MA0816.1', 'MA0691.1', 'MA0667.1', 'MA0665.1', 'MA0832.2', 'MA1619.2', 'MA1472.3', 'MA1641.2', 'MA0500.3', 'MA1100.3', 'MA0048.3', 'MA0698.2', 'MA1638.2', 'MA1467.3', 'MA1642.2', 'MA0668.3', 'MA1109.2', 'MA1123.3', 'MA1618.2', 'MA1593.2', 'MA1629.2', 'MA0697.3', 'MA1628.2', 'UN0112.1', 'UN0808.1', 'UN0331.2', 'UN0312.2', 'UN0676.2', 'UN0811.1', 'UN0616.2']
    motifs_cluster_006 = []
    for ma_id in cluster_006:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_006.append(motif)
    motifs_fetched.append(motifs_cluster_006)
    total_count=total_count+len(motifs_cluster_006)

    cluster_007 = ['MA1470.2', 'MA0794.1', 'MA0615.2', 'MA0862.1', 'MA0595.1', 'MA0596.1', 'MA0018.5', 'MA1951.2', 'MA0492.2', 'MA1133.2', 'MA1127.1', 'MA1140.3', 'MA0488.2', 'MA1632.2', 'MA0609.3', 'MA0656.2', 'MA0840.2', 'MA0605.3', 'MA1126.2', 'MA1131.2', 'MA1145.2', 'MA0834.2', 'MA1136.1', 'MA1129.1', 'MA1139.2', 'MA0604.1', 'MA1143.2', 'MA1475.2', 'MA1474.2', 'MA0839.2', 'MA1466.2', 'MA0638.2', 'MA0844.2', 'UN0110.2', 'UN0122.2', 'UN0655.2', 'UN0111.2']
    motifs_cluster_007 = []
    for ma_id in cluster_007:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_007.append(motif)
    motifs_fetched.append(motifs_cluster_007)
    total_count=total_count+len(motifs_cluster_007)

    cluster_008 = ['MA1603.2', 'MA0084.2', 'MA1152.2', 'MA0143.5', 'MA0514.3', 'MA1120.2', 'MA0867.3', 'MA0868.3', 'MA1563.2', 'MA0077.2', 'MA0087.3', 'MA1562.2', 'MA0515.1', 'MA0078.3', 'MA2095.1', 'MA0508.4', 'MA0442.3', 'MA0869.3', 'UN0814.1', 'UN0252.2', 'UN0262.2']
    motifs_cluster_008 = []
    for ma_id in cluster_008:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_008.append(motif)
    motifs_fetched.append(motifs_cluster_008)
    total_count=total_count+len(motifs_cluster_008)

    cluster_009 = ['MA1652.2', 'UN0606.2', 'UN0602.2', 'UN0804.1']
    motifs_cluster_009 = []
    for ma_id in cluster_009:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_009.append(motif)
    motifs_fetched.append(motifs_cluster_009)
    total_count=total_count+len(motifs_cluster_009)

    cluster_010 = ['MA0509.3', 'MA0798.3', 'MA0510.3', 'MA0600.3', 'MA0799.3', 'MA1554.2', 'MA1724.2']
    motifs_cluster_010 = []
    for ma_id in cluster_010:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_010.append(motif)
    motifs_fetched.append(motifs_cluster_010)
    total_count=total_count+len(motifs_cluster_010)

    cluster_011 = ['MA0602.2', 'MA1651.2', 'MA0095.4', 'UN0208.2', 'UN0241.2']
    motifs_cluster_011 = []
    for ma_id in cluster_011:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_011.append(motif)
    motifs_fetched.append(motifs_cluster_011)
    total_count=total_count+len(motifs_cluster_011)

    cluster_012 = ['MA0029.2', 'MA0766.3', 'MA0482.3', 'MA0037.5', 'MA1104.3', 'MA1970.2', 'MA0035.5', 'MA0036.4', 'UN0642.2']
    motifs_cluster_012 = []
    for ma_id in cluster_012:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_012.append(motif)
    motifs_fetched.append(motifs_cluster_012)
    total_count=total_count+len(motifs_cluster_012)

    cluster_013 = ['MA0145.2', 'MA1967.2', 'MA1945.2', 'MA1934.2', 'MA1941.2', 'MA1966.2', 'UN0509.1', 'UN0525.2', 'UN0507.2', 'UN0532.2', 'UN0529.2']
    motifs_cluster_013 = []
    for ma_id in cluster_013:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_013.append(motif)
    motifs_fetched.append(motifs_cluster_013)
    total_count=total_count+len(motifs_cluster_013)

    cluster_014 = ['MA2001.2', 'MA1118.2', 'MA1119.2', 'UN0231.2', 'UN0177.2']
    motifs_cluster_014 = []
    for ma_id in cluster_014:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_014.append(motif)
    motifs_fetched.append(motifs_cluster_014)
    total_count=total_count+len(motifs_cluster_014)

    cluster_015 = ['MA0142.1', 'MA1962.1', 'MA0787.1', 'MA0788.1', 'MA0507.3', 'MA0784.3', 'MA1115.2', 'MA0786.2', 'MA0627.3', 'MA0789.1', 'MA0785.2', 'MA0792.1', 'UN0577.1', 'UN0581.2']
    motifs_cluster_015 = []
    for ma_id in cluster_015:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_015.append(motif)
    motifs_fetched.append(motifs_cluster_015)
    total_count=total_count+len(motifs_cluster_015)

    cluster_016 = ['MA1116.2', 'MA0624.3', 'MA0625.3', 'MA0152.3', 'MA0606.3', 'MA1525.3', 'MA0130.1', 'UN0336.2', 'UN0175.2']
    motifs_cluster_016 = []
    for ma_id in cluster_016:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_016.append(motif)
    motifs_fetched.append(motifs_cluster_016)
    total_count=total_count+len(motifs_cluster_016)

    cluster_017 = ['MA0014.4', 'MA1995.2', 'MA0819.3', 'MA1108.3', 'MA0004.1', 'MA0058.4', 'MA0622.2', 'MA0825.2', 'MA1106.2', 'MA0006.2', 'MA0259.2', 'MA2325.1', 'MA0632.3', 'MA0829.3', 'MA0093.4', 'MA0663.1', 'MA0828.3', 'MA0620.4', 'MA0636.1', 'MA0526.5', 'MA0831.3', 'MA0147.4', 'MA0104.5', 'MA0626.2', 'MA0464.3', 'MA1560.2', 'MA1099.3', 'MA0664.2', 'MA0692.2', 'MA0871.3', 'MA0603.2', 'MA1464.2', 'MA0059.2', 'MA1493.1', 'MA0821.2', 'MA0823.1', 'MA0608.1', 'MA0616.3', 'MA0649.2', 'UN0580.2', 'UN0624.2', 'UN0803.1']
    motifs_cluster_017 = []
    for ma_id in cluster_017:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_017.append(motif)
    motifs_fetched.append(motifs_cluster_017)
    total_count=total_count+len(motifs_cluster_017)

    cluster_018 = ['MA1709.2', 'MA1418.2', 'MA0050.4', 'MA1623.2', 'MA0051.2', 'MA0517.2', 'MA0653.1', 'MA1420.1', 'MA0772.2', 'MA0652.2', 'MA1419.2', 'UN0152.2']
    motifs_cluster_018 = []
    for ma_id in cluster_018:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_018.append(motif)
    motifs_fetched.append(motifs_cluster_018)
    total_count=total_count+len(motifs_cluster_018)

    cluster_019 = ['MA0463.3', 'MA0731.1', 'MA0520.2', 'MA0518.2', 'MA0137.4', 'MA1624.2', 'MA0144.3', 'MA0519.2', 'MA1625.2', 'UN0663.2', 'UN0610.2']
    motifs_cluster_019 = []
    for ma_id in cluster_019:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_019.append(motif)
    motifs_fetched.append(motifs_cluster_019)
    total_count=total_count+len(motifs_cluster_019)

    cluster_020 = ['MA0139.2', 'MA0155.1', 'MA1102.3']
    motifs_cluster_020 = []
    for ma_id in cluster_020:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_020.append(motif)
    motifs_fetched.append(motifs_cluster_020)
    total_count=total_count+len(motifs_cluster_020)

    cluster_021 = ['MA1578.2', 'MA0528.3', 'MA1965.2', 'MA2328.1', 'MA0753.3', 'MA1627.2', 'MA1516.2', 'MA1107.3', 'MA0741.1', 'MA0746.3', 'MA1517.2', 'MA0747.2', 'MA1512.2', 'MA1564.2', 'MA0039.5', 'MA1513.2', 'MA1515.2', 'MA0493.3', 'MA1959.2', 'MA0599.1', 'MA0079.5', 'MA0685.2', 'MA0516.3', 'MA0740.2', 'MA0742.2', 'MA1511.2', 'MA1522.2', 'MA1961.2', 'MA1630.3', 'MA1653.2', 'UN0191.2', 'UN0628.1', 'UN0654.2', 'UN0800.1', 'UN0327.2', 'UN0326.2']
    motifs_cluster_021 = []
    for ma_id in cluster_021:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_021.append(motif)
    motifs_fetched.append(motifs_cluster_021)
    total_count=total_count+len(motifs_cluster_021)

    cluster_022 = ['MA1524.3', 'MA0461.3', 'MA0091.2', 'MA1570.1', 'MA0623.2', 'MA0669.1', 'MA0607.2', 'MA1568.2', 'MA1468.1', 'MA0817.2', 'MA0826.1', 'MA0818.2', 'MA0678.1', 'MA0827.1', 'UN0131.1', 'UN0143.2', 'UN0150.2']
    motifs_cluster_022 = []
    for ma_id in cluster_022:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_022.append(motif)
    motifs_fetched.append(motifs_cluster_022)
    total_count=total_count+len(motifs_cluster_022)

    cluster_023 = ['MA0748.3', 'MA1583.2', 'MA0865.3', 'MA0471.3', 'MA1122.2']
    motifs_cluster_023 = []
    for ma_id in cluster_023:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_023.append(motif)
    motifs_fetched.append(motifs_cluster_023)
    total_count=total_count+len(motifs_cluster_023)

    cluster_024 = ['MA0074.1', 'MA1533.2', 'MA2327.1', 'MA2323.1', 'MA0007.4', 'MA0113.4', 'MA0727.2']
    motifs_cluster_024 = []
    for ma_id in cluster_024:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_024.append(motif)
    motifs_fetched.append(motifs_cluster_024)
    total_count=total_count+len(motifs_cluster_024)

    cluster_025 = ['MA1657.2', 'MA0689.1', 'MA0806.1', 'MA0807.1', 'MA0801.1', 'MA0803.1', 'MA0805.1', 'MA0690.3', 'MA1567.3', 'MA0688.2', 'MA1565.2', 'MA0800.2', 'MA0802.2', 'MA1566.3', 'MA1960.2', 'UN0488.2', 'UN0566.1', 'UN0650.2', 'UN0317.2', 'UN0491.2', 'UN0564.2', 'UN0565.2', 'UN0485.1', 'UN0572.1', 'UN0573.1']
    motifs_cluster_025 = []
    for ma_id in cluster_025:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_025.append(motif)
    motifs_fetched.append(motifs_cluster_025)
    total_count=total_count+len(motifs_cluster_025)

    cluster_026 = ['MA1540.3', 'MA1541.2', 'MA0072.2', 'MA0071.1', 'MA1150.2', 'MA1151.2', 'MA1535.2', 'MA1536.2', 'MA2337.1', 'MA1110.3', 'MA1996.2', 'MA2324.1', 'MA0141.4', 'MA0505.3', 'MA0592.4', 'MA0643.2', 'MA0160.3', 'MA1112.3', 'MA1111.2', 'MA2338.1', 'UN0308.2']
    motifs_cluster_026 = []
    for ma_id in cluster_026:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_026.append(motif)
    motifs_fetched.append(motifs_cluster_026)
    total_count=total_count+len(motifs_cluster_026)

    cluster_027 = ['MA1647.3', 'MA1509.1', 'MA1971.2']
    motifs_cluster_027 = []
    for ma_id in cluster_027:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_027.append(motif)
    motifs_fetched.append(motifs_cluster_027)
    total_count=total_count+len(motifs_cluster_027)

    cluster_028 = ['MA2099.1', 'MA0003.5', 'MA0814.3', 'MA0812.2', 'MA1569.2', 'UN0174.1', 'UN0809.1', 'UN0160.2', 'UN0218.2']
    motifs_cluster_028 = []
    for ma_id in cluster_028:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_028.append(motif)
    motifs_fetched.append(motifs_cluster_028)
    total_count=total_count+len(motifs_cluster_028)

    cluster_029 = ['MA0506.3', 'MA1650.2', 'MA0162.5', 'MA0472.2', 'MA0732.2', 'MA0733.2']
    motifs_cluster_029 = []
    for ma_id in cluster_029:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_029.append(motif)
    motifs_fetched.append(motifs_cluster_029)
    total_count=total_count+len(motifs_cluster_029)

    cluster_030 = ['MA0066.2', 'MA0112.4', 'MA0258.2', 'UN0680.1']
    motifs_cluster_030 = []
    for ma_id in cluster_030:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_030.append(motif)
    motifs_fetched.append(motifs_cluster_030)
    total_count=total_count+len(motifs_cluster_030)

    cluster_031 = ['MA1546.2', 'MA0069.1', 'MA0779.2', 'MA0781.2', 'MA0067.3', 'MA2094.1', 'UN0562.2', 'UN0561.2', 'UN0563.1']
    motifs_cluster_031 = []
    for ma_id in cluster_031:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_031.append(motif)
    motifs_fetched.append(motifs_cluster_031)
    total_count=total_count+len(motifs_cluster_031)

    cluster_032 = ['MA1531.2', 'MA1532.2', 'MA0729.1', 'MA0857.1', 'MA0728.1', 'MA0859.2', 'MA0017.3', 'MA1148.2', 'MA0115.1', 'MA1494.2', 'MA0065.3', 'MA1574.2', 'MA0677.2', 'MA1537.2', 'MA0504.2', 'MA1550.2', 'MA0855.1', 'MA0512.2', 'MA0856.1']
    motifs_cluster_032 = []
    for ma_id in cluster_032:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_032.append(motif)
    motifs_fetched.append(motifs_cluster_032)
    total_count=total_count+len(motifs_cluster_032)

    cluster_033 = ['MA0883.2', 'MA1480.2', 'MA0719.2', 'MA0714.2', 'MA0467.3', 'MA0682.3', 'MA0711.2', 'MA0648.2', 'MA0891.2', 'MA0712.3', 'MA1547.2', 'UN0534.1', 'UN0134.2']
    motifs_cluster_033 = []
    for ma_id in cluster_033:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_033.append(motif)
    motifs_fetched.append(motifs_cluster_033)
    total_count=total_count+len(motifs_cluster_033)

    cluster_034 = ['MA1731.2', 'UN0238.2', 'UN0228.2', 'UN0677.1']
    motifs_cluster_034 = []
    for ma_id in cluster_034:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_034.append(motif)
    motifs_fetched.append(motifs_cluster_034)
    total_count=total_count+len(motifs_cluster_034)

    cluster_035 = ['MA1985.1', 'MA1555.1', 'MA1556.1', 'MA1552.2', 'MA1553.2']
    motifs_cluster_035 = []
    for ma_id in cluster_035:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_035.append(motif)
    motifs_fetched.append(motifs_cluster_035)
    total_count=total_count+len(motifs_cluster_035)

    cluster_036 = ['MA1471.2', 'MA0122.4', 'MA0724.1', 'MA0896.2', 'MA0898.2', 'MA1518.3', 'MA0853.2', 'MA0874.2', 'MA0895.2', 'MA1608.2', 'MA1963.2', 'MA1499.2', 'MA1504.2', 'MA1507.2', 'MA0723.3', 'MA0902.3', 'MA0675.2', 'MA0892.2', 'MA0132.3', 'MA0903.2', 'MA0707.3', 'MA0887.2', 'MA0888.2', 'MA0705.2', 'MA0900.3', 'MA0628.2', 'MA0704.2', 'MA0612.3', 'MA0886.2', 'MA1505.2', 'MA0904.3', 'MA1495.2', 'MA1498.3', 'MA0875.2', 'MA0876.2', 'MA0700.3', 'MA0027.3', 'MA1481.2', 'MA0630.2', 'MA0662.2', 'MA0634.2', 'MA0721.2', 'MA0716.2', 'MA0717.2', 'MA1577.2', 'MA0654.2', 'MA0720.2', 'MA0890.2', 'MA0718.2', 'MA0894.2', 'MA0125.2', 'MA0699.2', 'MA0882.2', 'MA0880.2', 'MA0881.2', 'MA0879.3', 'MA0708.3', 'MA0666.3', 'MA0709.2', 'MA0893.3', 'MA1496.2', 'MA1501.2', 'MA1497.2', 'MA1500.2', 'MA0910.3', 'MA1502.2', 'MA0674.2', 'MA0621.2', 'MA0702.3', 'MA0661.2', 'MA2119.1', 'MA0706.2', 'MA2093.1', 'MA0642.3', 'MA0618.2', 'MA1519.2', 'MA0889.2', 'MA0644.3', 'MA0075.4', 'MA0701.3', 'MA0710.2', 'MA0722.2', 'MA0725.2', 'MA0726.2', 'MA0885.3', 'MA1476.3', 'MA0158.2', 'MA0912.2', 'MA0703.3', 'MA0068.2', 'MA0658.2', 'MA0854.2', 'MA1463.2', 'MA0793.2', 'MA1530.2', 'MA1549.2']
    motifs_cluster_036 = []
    for ma_id in cluster_036:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_036.append(motif)
    motifs_fetched.append(motifs_cluster_036)
    total_count=total_count+len(motifs_cluster_036)

    cluster_037 = ['MA2334.1', 'UN0145.2', 'UN0638.2']
    motifs_cluster_037 = []
    for ma_id in cluster_037:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_037.append(motif)
    motifs_fetched.append(motifs_cluster_037)
    total_count=total_count+len(motifs_cluster_037)

    cluster_038 = ['MA0684.3', 'MA0511.2', 'MA0002.3', 'MA1989.2', 'UN0587.2']
    motifs_cluster_038 = []
    for ma_id in cluster_038:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_038.append(motif)
    motifs_fetched.append(motifs_cluster_038)
    total_count=total_count+len(motifs_cluster_038)

    cluster_039 = ['UN0592.2', 'UN0647.2']
    motifs_cluster_039 = []
    for ma_id in cluster_039:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_039.append(motif)
    motifs_fetched.append(motifs_cluster_039)
    total_count=total_count+len(motifs_cluster_039)

    cluster_040 = ['MA1943.2', 'MA1937.2', 'MA1958.2', 'MA1932.2', 'MA1948.2', 'MA1940.2', 'MA1944.2', 'MA1957.1', 'MA1931.1', 'MA1949.2', 'UN0496.2', 'UN0513.2', 'UN0559.2', 'UN0568.2', 'UN0531.1', 'UN0569.2']
    motifs_cluster_040 = []
    for ma_id in cluster_040:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_040.append(motif)
    motifs_fetched.append(motifs_cluster_040)
    total_count=total_count+len(motifs_cluster_040)

    cluster_041 = ['MA0502.3', 'MA2339.1', 'MA0060.4', 'MA1644.2', 'UN0819.1']
    motifs_cluster_041 = []
    for ma_id in cluster_041:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_041.append(motif)
    motifs_fetched.append(motifs_cluster_041)
    total_count=total_count+len(motifs_cluster_041)

    cluster_042 = ['MA0597.3', 'MA1726.2', 'UN0333.2', 'UN0630.2']
    motifs_cluster_042 = []
    for ma_id in cluster_042:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_042.append(motif)
    motifs_fetched.append(motifs_cluster_042)
    total_count=total_count+len(motifs_cluster_042)

    cluster_043 = ['MA0687.2', 'MA0080.7', 'MA0081.3', 'MA1508.2', 'MA1947.2', 'MA1935.2', 'MA1942.2', 'MA1952.2', 'MA1950.2', 'MA1954.2', 'MA1955.2', 'MA1953.2', 'MA1956.2', 'MA1936.2', 'MA1946.2', 'UN0508.1', 'UN0511.2']
    motifs_cluster_043 = []
    for ma_id in cluster_043:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_043.append(motif)
    motifs_fetched.append(motifs_cluster_043)
    total_count=total_count+len(motifs_cluster_043)

    cluster_044 = ['UN0597.1', 'UN0622.1', 'UN0667.2']
    motifs_cluster_044 = []
    for ma_id in cluster_044:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_044.append(motif)
    motifs_fetched.append(motifs_cluster_044)
    total_count=total_count+len(motifs_cluster_044)

    cluster_045 = ['MA1604.2', 'MA1637.2', 'MA0154.5', 'MA2122.1', 'MA0813.1', 'MA0815.1', 'MA0872.1', 'MA0810.2', 'MA0524.3', 'MA0811.2']
    motifs_cluster_045 = []
    for ma_id in cluster_045:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_045.append(motif)
    motifs_fetched.append(motifs_cluster_045)
    total_count=total_count+len(motifs_cluster_045)

    cluster_046 = ['MA0629.2', 'MA1730.2', 'MA1655.2', 'MA1725.2']
    motifs_cluster_046 = []
    for ma_id in cluster_046:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_046.append(motif)
    motifs_fetched.append(motifs_cluster_046)
    total_count=total_count+len(motifs_cluster_046)

    cluster_047 = ['MA0743.3', 'MA0744.3', 'MA1105.3', 'MA0647.2', 'MA1968.2']
    motifs_cluster_047 = []
    for ma_id in cluster_047:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_047.append(motif)
    motifs_fetched.append(motifs_cluster_047)
    total_count=total_count+len(motifs_cluster_047)

    cluster_048 = ['MA0111.1', 'MA0631.2']
    motifs_cluster_048 = []
    for ma_id in cluster_048:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_048.append(motif)
    motifs_fetched.append(motifs_cluster_048)
    total_count=total_count+len(motifs_cluster_048)

    cluster_049 = ['MA1579.2', 'MA1153.2', 'MA1964.2', 'UN0245.2']
    motifs_cluster_049 = []
    for ma_id in cluster_049:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_049.append(motif)
    motifs_fetched.append(motifs_cluster_049)
    total_count=total_count+len(motifs_cluster_049)

    cluster_050 = ['MA1649.2', 'MA2096.1']
    motifs_cluster_050 = []
    for ma_id in cluster_050:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_050.append(motif)
    motifs_fetched.append(motifs_cluster_050)
    total_count=total_count+len(motifs_cluster_050)

    cluster_051 = ['MA0135.2', 'MA0680.3', 'MA0780.1', 'MA0757.2', 'MA2102.1', 'MA0756.3', 'MA0679.3', 'MA0754.3', 'MA0755.2']
    motifs_cluster_051 = []
    for ma_id in cluster_051:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_051.append(motif)
    motifs_fetched.append(motifs_cluster_051)
    total_count=total_count+len(motifs_cluster_051)

    cluster_052 = ['MA0610.2', 'MA1478.2', 'MA1479.2', 'MA1707.2']
    motifs_cluster_052 = []
    for ma_id in cluster_052:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_052.append(motif)
    motifs_fetched.append(motifs_cluster_052)
    total_count=total_count+len(motifs_cluster_052)

    cluster_053 = ['UN0494.1', 'UN0502.1', 'UN0500.1', 'UN0501.1']
    motifs_cluster_053 = []
    for ma_id in cluster_053:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_053.append(motif)
    motifs_fetched.append(motifs_cluster_053)
    total_count=total_count+len(motifs_cluster_053)

    cluster_054 = ['MA0146.3', 'UN0658.2']
    motifs_cluster_054 = []
    for ma_id in cluster_054:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_054.append(motif)
    motifs_fetched.append(motifs_cluster_054)
    total_count=total_count+len(motifs_cluster_054)

    cluster_055 = ['UN0159.2', 'UN0593.2', 'UN0668.2']
    motifs_cluster_055 = []
    for ma_id in cluster_055:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_055.append(motif)
    motifs_fetched.append(motifs_cluster_055)
    total_count=total_count+len(motifs_cluster_055)

    cluster_056 = ['MA0527.2', 'MA0749.2', 'UN0185.2']
    motifs_cluster_056 = []
    for ma_id in cluster_056:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_056.append(motif)
    motifs_fetched.append(motifs_cluster_056)
    total_count=total_count+len(motifs_cluster_056)

    cluster_057 = ['MA0641.1', 'MA0761.3', 'MA2332.1', 'MA0062.4', 'MA0474.4', 'MA0598.4', 'MA0640.3', 'MA0473.4', 'MA1992.2', 'MA0645.2', 'MA1708.2', 'MA2329.1', 'MA0076.3', 'MA0750.3', 'MA2340.1', 'MA0762.2', 'MA1484.2', 'MA0098.4', 'MA0475.3', 'MA0763.2', 'MA0156.4', 'MA0760.2', 'MA0764.4', 'MA0765.4', 'MA0028.3', 'MA0759.3', 'MA0686.2', 'MA1483.3', 'UN0530.2', 'UN0536.2', 'UN0514.2', 'UN0516.2', 'UN0325.2', 'UN0340.2', 'UN0801.1']
    motifs_cluster_057 = []
    for ma_id in cluster_057:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_057.append(motif)
    motifs_fetched.append(motifs_cluster_057)
    total_count=total_count+len(motifs_cluster_057)

    cluster_058 = ['UN0584.2', 'UN0585.2', 'UN0583.2', 'UN0582.2', 'UN0586.2']
    motifs_cluster_058 = []
    for ma_id in cluster_058:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_058.append(motif)
    motifs_fetched.append(motifs_cluster_058)
    total_count=total_count+len(motifs_cluster_058)

    cluster_059 = ['MA0149.1', 'UN0596.2']
    motifs_cluster_059 = []
    for ma_id in cluster_059:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_059.append(motif)
    motifs_fetched.append(motifs_cluster_059)
    total_count=total_count+len(motifs_cluster_059)

    cluster_060 = ['MA0136.4', 'MA2326.1', 'UN0649.2']
    motifs_cluster_060 = []
    for ma_id in cluster_060:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_060.append(motif)
    motifs_fetched.append(motifs_cluster_060)
    total_count=total_count+len(motifs_cluster_060)

    cluster_061 = ['MA0131.3', 'UN0627.2']
    motifs_cluster_061 = []
    for ma_id in cluster_061:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_061.append(motif)
    motifs_fetched.append(motifs_cluster_061)
    total_count=total_count+len(motifs_cluster_061)

    cluster_062 = ['MA0108.3', 'MA0151.1', 'MA0601.2']
    motifs_cluster_062 = []
    for ma_id in cluster_062:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_062.append(motif)
    motifs_fetched.append(motifs_cluster_062)
    total_count=total_count+len(motifs_cluster_062)

    cluster_063 = ['UN0515.2', 'UN0535.2', 'UN0541.2', 'UN0552.2', 'UN0543.1', 'UN0544.2']
    motifs_cluster_063 = []
    for ma_id in cluster_063:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_063.append(motif)
    motifs_fetched.append(motifs_cluster_063)
    total_count=total_count+len(motifs_cluster_063)

    cluster_064 = ['MA0758.1', 'MA0469.4', 'MA0024.3', 'MA0470.3', 'MA0864.3']
    motifs_cluster_064 = []
    for ma_id in cluster_064:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_064.append(motif)
    motifs_fetched.append(motifs_cluster_064)
    total_count=total_count+len(motifs_cluster_064)

    cluster_065 = ['MA2002.2', 'UN0227.2']
    motifs_cluster_065 = []
    for ma_id in cluster_065:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_065.append(motif)
    motifs_fetched.append(motifs_cluster_065)
    total_count=total_count+len(motifs_cluster_065)

    cluster_066 = ['MA0101.1', 'MA0107.1', 'MA0105.4', 'MA0778.2']
    motifs_cluster_066 = []
    for ma_id in cluster_066:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_066.append(motif)
    motifs_fetched.append(motifs_cluster_066)
    total_count=total_count+len(motifs_cluster_066)

    cluster_067 = ['UN0305.2', 'UN0645.2']
    motifs_cluster_067 = []
    for ma_id in cluster_067:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_067.append(motif)
    motifs_fetched.append(motifs_cluster_067)
    total_count=total_count+len(motifs_cluster_067)

    cluster_068 = ['UN0518.2', 'UN0560.2', 'UN0567.2']
    motifs_cluster_068 = []
    for ma_id in cluster_068:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_068.append(motif)
    motifs_fetched.append(motifs_cluster_068)
    total_count=total_count+len(motifs_cluster_068)

    cluster_069 = ['MA1645.2', 'MA0063.3', 'MA1523.2', 'MA1994.2', 'MA0672.2', 'MA0673.2', 'MA2003.2', 'MA0124.3', 'MA0914.2']
    motifs_cluster_069 = []
    for ma_id in cluster_069:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_069.append(motif)
    motifs_fetched.append(motifs_cluster_069)
    total_count=total_count+len(motifs_cluster_069)

    cluster_070 = ['MA0861.2', 'MA0106.3', 'MA0525.2', 'UN0130.1']
    motifs_cluster_070 = []
    for ma_id in cluster_070:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_070.append(motif)
    motifs_fetched.append(motifs_cluster_070)
    total_count=total_count+len(motifs_cluster_070)

    cluster_071 = ['UN0171.1', 'UN0316.1']
    motifs_cluster_071 = []
    for ma_id in cluster_071:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_071.append(motif)
    motifs_fetched.append(motifs_cluster_071)
    total_count=total_count+len(motifs_cluster_071)

    cluster_072 = ['MA2341.1', 'UN0653.2', 'UN0659.2', 'UN0605.2']
    motifs_cluster_072 = []
    for ma_id in cluster_072:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_072.append(motif)
    motifs_fetched.append(motifs_cluster_072)
    total_count=total_count+len(motifs_cluster_072)

    cluster_073 = ['MA1938.2', 'UN0542.1', 'UN0517.2']
    motifs_cluster_073 = []
    for ma_id in cluster_073:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_073.append(motif)
    motifs_fetched.append(motifs_cluster_073)
    total_count=total_count+len(motifs_cluster_073)

    cluster_074 = ['MA0046.3', 'MA0153.2', 'MA0683.2', 'MA0790.2', 'MA0791.2']
    motifs_cluster_074 = []
    for ma_id in cluster_074:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_074.append(motif)
    motifs_fetched.append(motifs_cluster_074)
    total_count=total_count+len(motifs_cluster_074)

    cluster_075 = ['UN0324.2', 'UN0820.1']
    motifs_cluster_075 = []
    for ma_id in cluster_075:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_075.append(motif)
    motifs_fetched.append(motifs_cluster_075)
    total_count=total_count+len(motifs_cluster_075)

    cluster_076 = ['MA0842.3', 'MA0117.3', 'MA0495.4', 'UN0197.2', 'UN0319.2']
    motifs_cluster_076 = []
    for ma_id in cluster_076:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_076.append(motif)
    motifs_fetched.append(motifs_cluster_076)
    total_count=total_count+len(motifs_cluster_076)

    cluster_077 = ['MA1727.2', 'MA2126.1', 'UN0639.2']
    motifs_cluster_077 = []
    for ma_id in cluster_077:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_077.append(motif)
    motifs_fetched.append(motifs_cluster_077)
    total_count=total_count+len(motifs_cluster_077)

    cluster_078 = ['MA1587.1', 'MA1596.1']
    motifs_cluster_078 = []
    for ma_id in cluster_078:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_078.append(motif)
    motifs_fetched.append(motifs_cluster_078)
    total_count=total_count+len(motifs_cluster_078)

    cluster_079 = ['MA1584.2', 'MA0696.1', 'MA0751.2', 'MA1491.3', 'MA0735.2', 'MA0736.1', 'MA0737.1']
    motifs_cluster_079 = []
    for ma_id in cluster_079:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_079.append(motif)
    motifs_fetched.append(motifs_cluster_079)
    total_count=total_count+len(motifs_cluster_079)

    cluster_080 = ['UN0330.1', 'UN0812.1']
    motifs_cluster_080 = []
    for ma_id in cluster_080:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_080.append(motif)
    motifs_fetched.append(motifs_cluster_080)
    total_count=total_count+len(motifs_cluster_080)

    cluster_081 = ['MA0019.2', 'UN0611.2']
    motifs_cluster_081 = []
    for ma_id in cluster_081:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_081.append(motif)
    motifs_fetched.append(motifs_cluster_081)
    total_count=total_count+len(motifs_cluster_081)

    cluster_082 = ['UN0198.2', 'UN0232.2']
    motifs_cluster_082 = []
    for ma_id in cluster_082:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_082.append(motif)
    motifs_fetched.append(motifs_cluster_082)
    total_count=total_count+len(motifs_cluster_082)

    cluster_083 = ['MA0676.1', 'MA0114.5', 'MA0484.3']
    motifs_cluster_083 = []
    for ma_id in cluster_083:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_083.append(motif)
    motifs_fetched.append(motifs_cluster_083)
    total_count=total_count+len(motifs_cluster_083)

    cluster_084 = ['MA2330.1', 'UN0588.1']
    motifs_cluster_084 = []
    for ma_id in cluster_084:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_084.append(motif)
    motifs_fetched.append(motifs_cluster_084)
    total_count=total_count+len(motifs_cluster_084)

    cluster_085 = ['MA2098.1', 'UN0634.2']
    motifs_cluster_085 = []
    for ma_id in cluster_085:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_085.append(motif)
    motifs_fetched.append(motifs_cluster_085)
    total_count=total_count+len(motifs_cluster_085)

    cluster_086 = ['MA0768.3', 'MA0769.3', 'MA1421.1', 'MA0523.2', 'MA1991.2']
    motifs_cluster_086 = []
    for ma_id in cluster_086:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_086.append(motif)
    motifs_fetched.append(motifs_cluster_086)
    total_count=total_count+len(motifs_cluster_086)

    cluster_087 = ['UN0665.2', 'UN0681.2']
    motifs_cluster_087 = []
    for ma_id in cluster_087:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_087.append(motif)
    motifs_fetched.append(motifs_cluster_087)
    total_count=total_count+len(motifs_cluster_087)

    cluster_088 = ['MA0681.3', 'MA0713.1', 'MA0715.1', 'MA0611.3', 'MA0884.2', 'MA0468.1', 'UN0671.2']
    motifs_cluster_088 = []
    for ma_id in cluster_088:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_088.append(motif)
    motifs_fetched.append(motifs_cluster_088)
    total_count=total_count+len(motifs_cluster_088)

    cluster_089 = ['MA0070.2', 'MA1639.2', 'MA1113.3', 'MA1640.2']
    motifs_cluster_089 = []
    for ma_id in cluster_089:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_089.append(motif)
    motifs_fetched.append(motifs_cluster_089)
    total_count=total_count+len(motifs_cluster_089)

    cluster_090 = ['UN0651.2', 'UN0684.2']
    motifs_cluster_090 = []
    for ma_id in cluster_090:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_090.append(motif)
    motifs_fetched.append(motifs_cluster_090)
    total_count=total_count+len(motifs_cluster_090)

    cluster_091 = ['MA1602.2', 'MA0795.1', 'MA1557.1']
    motifs_cluster_091 = []
    for ma_id in cluster_091:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_091.append(motif)
    motifs_fetched.append(motifs_cluster_091)
    total_count=total_count+len(motifs_cluster_091)

    cluster_092 = ['MA0695.2', 'MA0694.2', 'MA0734.4', 'MA1990.2']
    motifs_cluster_092 = []
    for ma_id in cluster_092:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_092.append(motif)
    motifs_fetched.append(motifs_cluster_092)
    total_count=total_count+len(motifs_cluster_092)

    cluster_093 = ['MA1601.2', 'MA2097.1']
    motifs_cluster_093 = []
    for ma_id in cluster_093:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_093.append(motif)
    motifs_fetched.append(motifs_cluster_093)
    total_count=total_count+len(motifs_cluster_093)

    cluster_094 = ['MA1599.2', 'UN0662.2']
    motifs_cluster_094 = []
    for ma_id in cluster_094:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_094.append(motif)
    motifs_fetched.append(motifs_cluster_094)
    total_count=total_count+len(motifs_cluster_094)

    cluster_095 = ['MA0494.2', 'MA1969.2', 'MA0860.1', 'MA1576.2']
    motifs_cluster_095 = []
    for ma_id in cluster_095:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_095.append(motif)
    motifs_fetched.append(motifs_cluster_095)
    total_count=total_count+len(motifs_cluster_095)

    cluster_096 = ['MA1121.2', 'MA0090.4', 'MA0808.1', 'MA0809.3', 'UN0617.2']
    motifs_cluster_096 = []
    for ma_id in cluster_096:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_096.append(motif)
    motifs_fetched.append(motifs_cluster_096)
    total_count=total_count+len(motifs_cluster_096)

    cluster_097 = ['MA0056.3', 'MA0479.2']
    motifs_cluster_097 = []
    for ma_id in cluster_097:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_097.append(motif)
    motifs_fetched.append(motifs_cluster_097)
    total_count=total_count+len(motifs_cluster_097)

    cluster_098 = ['UN0524.2', 'UN0506.2', 'UN0512.2', 'UN0495.2', 'UN0503.2', 'UN0520.2', 'UN0528.2']
    motifs_cluster_098 = []
    for ma_id in cluster_098:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_098.append(motif)
    motifs_fetched.append(motifs_cluster_098)
    total_count=total_count+len(motifs_cluster_098)

    cluster_099 = ['UN0155.1', 'UN0631.1']
    motifs_cluster_099 = []
    for ma_id in cluster_099:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_099.append(motif)
    motifs_fetched.append(motifs_cluster_099)
    total_count=total_count+len(motifs_cluster_099)

    cluster_100 = ['MA1548.2', 'MA1615.2']
    motifs_cluster_100 = []
    for ma_id in cluster_100:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_100.append(motif)
    motifs_fetched.append(motifs_cluster_100)
    total_count=total_count+len(motifs_cluster_100)

    cluster_101 = ['MA0660.1', 'MA0773.1', 'MA0052.5', 'MA0497.2']
    motifs_cluster_101 = []
    for ma_id in cluster_101:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_101.append(motif)
    motifs_fetched.append(motifs_cluster_101)
    total_count=total_count+len(motifs_cluster_101)

    cluster_102 = ['MA0088.2', 'MA1716.2']
    motifs_cluster_102 = []
    for ma_id in cluster_102:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_102.append(motif)
    motifs_fetched.append(motifs_cluster_102)
    total_count=total_count+len(motifs_cluster_102)

    cluster_103 = ['MA1538.1', 'MA1539.1', 'MA1146.2', 'MA1147.2']
    motifs_cluster_103 = []
    for ma_id in cluster_103:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_103.append(motif)
    motifs_fetched.append(motifs_cluster_103)
    total_count=total_count+len(motifs_cluster_103)

    cluster_104 = ['MA0693.4', 'MA1534.2']
    motifs_cluster_104 = []
    for ma_id in cluster_104:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_104.append(motif)
    motifs_fetched.append(motifs_cluster_104)
    total_count=total_count+len(motifs_cluster_104)

    cluster_105 = ['MA1542.2', 'MA1646.2']
    motifs_cluster_105 = []
    for ma_id in cluster_105:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_105.append(motif)
    motifs_fetched.append(motifs_cluster_105)
    total_count=total_count+len(motifs_cluster_105)

    cluster_106 = ['MA0038.3', 'MA0483.2']
    motifs_cluster_106 = []
    for ma_id in cluster_106:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_106.append(motif)
    motifs_fetched.append(motifs_cluster_106)
    total_count=total_count+len(motifs_cluster_106)

    cluster_107 = ['MA1544.2', 'UN0806.1']
    motifs_cluster_107 = []
    for ma_id in cluster_107:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_107.append(motif)
    motifs_fetched.append(motifs_cluster_107)
    total_count=total_count+len(motifs_cluster_107)

    cluster_108 = ['UN0578.1', 'UN0579.2']
    motifs_cluster_108 = []
    for ma_id in cluster_108:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_108.append(motif)
    motifs_fetched.append(motifs_cluster_108)
    total_count=total_count+len(motifs_cluster_108)

    cluster_109 = ['MA0866.1', 'MA0870.1']
    motifs_cluster_109 = []
    for ma_id in cluster_109:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_109.append(motif)
    motifs_fetched.append(motifs_cluster_109)
    total_count=total_count+len(motifs_cluster_109)

    cluster_110 = ['MA0730.1', 'MA0858.1', 'MA0159.1', 'MA1149.2']
    motifs_cluster_110 = []
    for ma_id in cluster_110:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_110.append(motif)
    motifs_fetched.append(motifs_cluster_110)
    total_count=total_count+len(motifs_cluster_110)

    cluster_111 = ['MA1155.1', 'UN0243.2']
    motifs_cluster_111 = []
    for ma_id in cluster_111:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_111.append(motif)
    motifs_fetched.append(motifs_cluster_111)
    total_count=total_count+len(motifs_cluster_111)

    cluster_112 = ['MA0009.2', 'MA0804.2']
    motifs_cluster_112 = []
    for ma_id in cluster_112:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_112.append(motif)
    motifs_fetched.append(motifs_cluster_112)
    total_count=total_count+len(motifs_cluster_112)

    cluster_113 = ['UN0550.2', 'UN0555.2']
    motifs_cluster_113 = []
    for ma_id in cluster_113:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_113.append(motif)
    motifs_fetched.append(motifs_cluster_113)
    total_count=total_count+len(motifs_cluster_113)

    cluster_114 = ['MA0771.1', 'MA0486.2', 'MA0770.1']
    motifs_cluster_114 = []
    for ma_id in cluster_114:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_114.append(motif)
    motifs_fetched.append(motifs_cluster_114)
    total_count=total_count+len(motifs_cluster_114)

    cluster_115 = ['MA0073.2']
    motifs_cluster_115 = []
    for ma_id in cluster_115:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_115.append(motif)
    motifs_fetched.append(motifs_cluster_115)
    total_count=total_count+len(motifs_cluster_115)

    cluster_116 = ['MA0083.3']
    motifs_cluster_116 = []
    for ma_id in cluster_116:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_116.append(motif)
    motifs_fetched.append(motifs_cluster_116)
    total_count=total_count+len(motifs_cluster_116)

    cluster_117 = ['MA0092.2']
    motifs_cluster_117 = []
    for ma_id in cluster_117:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_117.append(motif)
    motifs_fetched.append(motifs_cluster_117)
    total_count=total_count+len(motifs_cluster_117)

    cluster_118 = ['MA0116.1']
    motifs_cluster_118 = []
    for ma_id in cluster_118:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_118.append(motif)
    motifs_fetched.append(motifs_cluster_118)
    total_count=total_count+len(motifs_cluster_118)

    cluster_119 = ['MA0138.3']
    motifs_cluster_119 = []
    for ma_id in cluster_119:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_119.append(motif)
    motifs_fetched.append(motifs_cluster_119)
    total_count=total_count+len(motifs_cluster_119)

    cluster_120 = ['MA0140.3']
    motifs_cluster_120 = []
    for ma_id in cluster_120:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_120.append(motif)
    motifs_fetched.append(motifs_cluster_120)
    total_count=total_count+len(motifs_cluster_120)

    cluster_121 = ['MA0163.1']
    motifs_cluster_121 = []
    for ma_id in cluster_121:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_121.append(motif)
    motifs_fetched.append(motifs_cluster_121)
    total_count=total_count+len(motifs_cluster_121)

    cluster_122 = ['MA0164.2']
    motifs_cluster_122 = []
    for ma_id in cluster_122:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_122.append(motif)
    motifs_fetched.append(motifs_cluster_122)
    total_count=total_count+len(motifs_cluster_122)

    cluster_123 = ['MA0594.3']
    motifs_cluster_123 = []
    for ma_id in cluster_123:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_123.append(motif)
    motifs_fetched.append(motifs_cluster_123)
    total_count=total_count+len(motifs_cluster_123)

    cluster_124 = ['MA0619.2']
    motifs_cluster_124 = []
    for ma_id in cluster_124:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_124.append(motif)
    motifs_fetched.append(motifs_cluster_124)
    total_count=total_count+len(motifs_cluster_124)

    cluster_125 = ['MA0657.2']
    motifs_cluster_125 = []
    for ma_id in cluster_125:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_125.append(motif)
    motifs_fetched.append(motifs_cluster_125)
    total_count=total_count+len(motifs_cluster_125)

    cluster_126 = ['MA0752.2']
    motifs_cluster_126 = []
    for ma_id in cluster_126:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_126.append(motif)
    motifs_fetched.append(motifs_cluster_126)
    total_count=total_count+len(motifs_cluster_126)

    cluster_127 = ['MA0776.1']
    motifs_cluster_127 = []
    for ma_id in cluster_127:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_127.append(motif)
    motifs_fetched.append(motifs_cluster_127)
    total_count=total_count+len(motifs_cluster_127)

    cluster_128 = ['MA0777.1']
    motifs_cluster_128 = []
    for ma_id in cluster_128:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_128.append(motif)
    motifs_fetched.append(motifs_cluster_128)
    total_count=total_count+len(motifs_cluster_128)

    cluster_129 = ['MA0863.1']
    motifs_cluster_129 = []
    for ma_id in cluster_129:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_129.append(motif)
    motifs_fetched.append(motifs_cluster_129)
    total_count=total_count+len(motifs_cluster_129)

    cluster_130 = ['MA0897.2']
    motifs_cluster_130 = []
    for ma_id in cluster_130:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_130.append(motif)
    motifs_fetched.append(motifs_cluster_130)
    total_count=total_count+len(motifs_cluster_130)

    cluster_131 = ['MA1117.2']
    motifs_cluster_131 = []
    for ma_id in cluster_131:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_131.append(motif)
    motifs_fetched.append(motifs_cluster_131)
    total_count=total_count+len(motifs_cluster_131)

    cluster_132 = ['MA1124.1']
    motifs_cluster_132 = []
    for ma_id in cluster_132:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_132.append(motif)
    motifs_fetched.append(motifs_cluster_132)
    total_count=total_count+len(motifs_cluster_132)

    cluster_133 = ['MA1125.2']
    motifs_cluster_133 = []
    for ma_id in cluster_133:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_133.append(motif)
    motifs_fetched.append(motifs_cluster_133)
    total_count=total_count+len(motifs_cluster_133)

    cluster_134 = ['MA1154.2']
    motifs_cluster_134 = []
    for ma_id in cluster_134:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_134.append(motif)
    motifs_fetched.append(motifs_cluster_134)
    total_count=total_count+len(motifs_cluster_134)

    cluster_135 = ['MA1514.2']
    motifs_cluster_135 = []
    for ma_id in cluster_135:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_135.append(motif)
    motifs_fetched.append(motifs_cluster_135)
    total_count=total_count+len(motifs_cluster_135)

    cluster_136 = ['MA1529.2']
    motifs_cluster_136 = []
    for ma_id in cluster_136:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_136.append(motif)
    motifs_fetched.append(motifs_cluster_136)
    total_count=total_count+len(motifs_cluster_136)

    cluster_137 = ['MA1545.2']
    motifs_cluster_137 = []
    for ma_id in cluster_137:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_137.append(motif)
    motifs_fetched.append(motifs_cluster_137)
    total_count=total_count+len(motifs_cluster_137)

    cluster_138 = ['MA1561.2']
    motifs_cluster_138 = []
    for ma_id in cluster_138:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_138.append(motif)
    motifs_fetched.append(motifs_cluster_138)
    total_count=total_count+len(motifs_cluster_138)

    cluster_139 = ['MA1573.2']
    motifs_cluster_139 = []
    for ma_id in cluster_139:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_139.append(motif)
    motifs_fetched.append(motifs_cluster_139)
    total_count=total_count+len(motifs_cluster_139)

    cluster_140 = ['MA1575.2']
    motifs_cluster_140 = []
    for ma_id in cluster_140:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_140.append(motif)
    motifs_fetched.append(motifs_cluster_140)
    total_count=total_count+len(motifs_cluster_140)

    cluster_141 = ['MA1580.1']
    motifs_cluster_141 = []
    for ma_id in cluster_141:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_141.append(motif)
    motifs_fetched.append(motifs_cluster_141)
    total_count=total_count+len(motifs_cluster_141)

    cluster_142 = ['MA1581.2']
    motifs_cluster_142 = []
    for ma_id in cluster_142:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_142.append(motif)
    motifs_fetched.append(motifs_cluster_142)
    total_count=total_count+len(motifs_cluster_142)

    cluster_143 = ['MA1585.2']
    motifs_cluster_143 = []
    for ma_id in cluster_143:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_143.append(motif)
    motifs_fetched.append(motifs_cluster_143)
    total_count=total_count+len(motifs_cluster_143)

    cluster_144 = ['MA1588.1']
    motifs_cluster_144 = []
    for ma_id in cluster_144:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_144.append(motif)
    motifs_fetched.append(motifs_cluster_144)
    total_count=total_count+len(motifs_cluster_144)

    cluster_145 = ['MA1589.2']
    motifs_cluster_145 = []
    for ma_id in cluster_145:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_145.append(motif)
    motifs_fetched.append(motifs_cluster_145)
    total_count=total_count+len(motifs_cluster_145)

    cluster_146 = ['MA1594.1']
    motifs_cluster_146 = []
    for ma_id in cluster_146:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_146.append(motif)
    motifs_fetched.append(motifs_cluster_146)
    total_count=total_count+len(motifs_cluster_146)

    cluster_147 = ['MA1597.1']
    motifs_cluster_147 = []
    for ma_id in cluster_147:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_147.append(motif)
    motifs_fetched.append(motifs_cluster_147)
    total_count=total_count+len(motifs_cluster_147)

    cluster_148 = ['MA1600.2']
    motifs_cluster_148 = []
    for ma_id in cluster_148:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_148.append(motif)
    motifs_fetched.append(motifs_cluster_148)
    total_count=total_count+len(motifs_cluster_148)

    cluster_149 = ['MA1616.2']
    motifs_cluster_149 = []
    for ma_id in cluster_149:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_149.append(motif)
    motifs_fetched.append(motifs_cluster_149)
    total_count=total_count+len(motifs_cluster_149)

    cluster_150 = ['MA1654.2']
    motifs_cluster_150 = []
    for ma_id in cluster_150:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_150.append(motif)
    motifs_fetched.append(motifs_cluster_150)
    total_count=total_count+len(motifs_cluster_150)

    cluster_151 = ['MA1656.2']
    motifs_cluster_151 = []
    for ma_id in cluster_151:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_151.append(motif)
    motifs_fetched.append(motifs_cluster_151)
    total_count=total_count+len(motifs_cluster_151)

    cluster_152 = ['MA1684.1']
    motifs_cluster_152 = []
    for ma_id in cluster_152:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_152.append(motif)
    motifs_fetched.append(motifs_cluster_152)
    total_count=total_count+len(motifs_cluster_152)

    cluster_153 = ['MA1710.2']
    motifs_cluster_153 = []
    for ma_id in cluster_153:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_153.append(motif)
    motifs_fetched.append(motifs_cluster_153)
    total_count=total_count+len(motifs_cluster_153)

    cluster_154 = ['MA1711.2']
    motifs_cluster_154 = []
    for ma_id in cluster_154:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_154.append(motif)
    motifs_fetched.append(motifs_cluster_154)
    total_count=total_count+len(motifs_cluster_154)

    cluster_155 = ['MA1712.2']
    motifs_cluster_155 = []
    for ma_id in cluster_155:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_155.append(motif)
    motifs_fetched.append(motifs_cluster_155)
    total_count=total_count+len(motifs_cluster_155)

    cluster_156 = ['MA1713.2']
    motifs_cluster_156 = []
    for ma_id in cluster_156:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_156.append(motif)
    motifs_fetched.append(motifs_cluster_156)
    total_count=total_count+len(motifs_cluster_156)

    cluster_157 = ['MA1714.2']
    motifs_cluster_157 = []
    for ma_id in cluster_157:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_157.append(motif)
    motifs_fetched.append(motifs_cluster_157)
    total_count=total_count+len(motifs_cluster_157)

    cluster_158 = ['MA1715.1']
    motifs_cluster_158 = []
    for ma_id in cluster_158:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_158.append(motif)
    motifs_fetched.append(motifs_cluster_158)
    total_count=total_count+len(motifs_cluster_158)

    cluster_159 = ['MA1717.2']
    motifs_cluster_159 = []
    for ma_id in cluster_159:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_159.append(motif)
    motifs_fetched.append(motifs_cluster_159)
    total_count=total_count+len(motifs_cluster_159)

    cluster_160 = ['MA1718.1']
    motifs_cluster_160 = []
    for ma_id in cluster_160:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_160.append(motif)
    motifs_fetched.append(motifs_cluster_160)
    total_count=total_count+len(motifs_cluster_160)

    cluster_161 = ['MA1719.2']
    motifs_cluster_161 = []
    for ma_id in cluster_161:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_161.append(motif)
    motifs_fetched.append(motifs_cluster_161)
    total_count=total_count+len(motifs_cluster_161)

    cluster_162 = ['MA1720.2']
    motifs_cluster_162 = []
    for ma_id in cluster_162:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_162.append(motif)
    motifs_fetched.append(motifs_cluster_162)
    total_count=total_count+len(motifs_cluster_162)

    cluster_163 = ['MA1721.2']
    motifs_cluster_163 = []
    for ma_id in cluster_163:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_163.append(motif)
    motifs_fetched.append(motifs_cluster_163)
    total_count=total_count+len(motifs_cluster_163)

    cluster_164 = ['MA1722.2']
    motifs_cluster_164 = []
    for ma_id in cluster_164:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_164.append(motif)
    motifs_fetched.append(motifs_cluster_164)
    total_count=total_count+len(motifs_cluster_164)

    cluster_165 = ['MA1723.2']
    motifs_cluster_165 = []
    for ma_id in cluster_165:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_165.append(motif)
    motifs_fetched.append(motifs_cluster_165)
    total_count=total_count+len(motifs_cluster_165)

    cluster_166 = ['MA1728.2']
    motifs_cluster_166 = []
    for ma_id in cluster_166:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_166.append(motif)
    motifs_fetched.append(motifs_cluster_166)
    total_count=total_count+len(motifs_cluster_166)

    cluster_167 = ['MA1729.2']
    motifs_cluster_167 = []
    for ma_id in cluster_167:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_167.append(motif)
    motifs_fetched.append(motifs_cluster_167)
    total_count=total_count+len(motifs_cluster_167)

    cluster_168 = ['MA1929.2']
    motifs_cluster_168 = []
    for ma_id in cluster_168:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_168.append(motif)
    motifs_fetched.append(motifs_cluster_168)
    total_count=total_count+len(motifs_cluster_168)

    cluster_169 = ['MA1930.2']
    motifs_cluster_169 = []
    for ma_id in cluster_169:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_169.append(motif)
    motifs_fetched.append(motifs_cluster_169)
    total_count=total_count+len(motifs_cluster_169)

    cluster_170 = ['MA1972.1']
    motifs_cluster_170 = []
    for ma_id in cluster_170:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_170.append(motif)
    motifs_fetched.append(motifs_cluster_170)
    total_count=total_count+len(motifs_cluster_170)

    cluster_171 = ['MA1973.2']
    motifs_cluster_171 = []
    for ma_id in cluster_171:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_171.append(motif)
    motifs_fetched.append(motifs_cluster_171)
    total_count=total_count+len(motifs_cluster_171)

    cluster_172 = ['MA1974.2']
    motifs_cluster_172 = []
    for ma_id in cluster_172:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_172.append(motif)
    motifs_fetched.append(motifs_cluster_172)
    total_count=total_count+len(motifs_cluster_172)

    cluster_173 = ['MA1975.2']
    motifs_cluster_173 = []
    for ma_id in cluster_173:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_173.append(motif)
    motifs_fetched.append(motifs_cluster_173)
    total_count=total_count+len(motifs_cluster_173)

    cluster_174 = ['MA1976.2']
    motifs_cluster_174 = []
    for ma_id in cluster_174:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_174.append(motif)
    motifs_fetched.append(motifs_cluster_174)
    total_count=total_count+len(motifs_cluster_174)

    cluster_175 = ['MA1977.2']
    motifs_cluster_175 = []
    for ma_id in cluster_175:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_175.append(motif)
    motifs_fetched.append(motifs_cluster_175)
    total_count=total_count+len(motifs_cluster_175)

    cluster_176 = ['MA1978.2']
    motifs_cluster_176 = []
    for ma_id in cluster_176:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_176.append(motif)
    motifs_fetched.append(motifs_cluster_176)
    total_count=total_count+len(motifs_cluster_176)

    cluster_177 = ['MA1979.2']
    motifs_cluster_177 = []
    for ma_id in cluster_177:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_177.append(motif)
    motifs_fetched.append(motifs_cluster_177)
    total_count=total_count+len(motifs_cluster_177)

    cluster_178 = ['MA1980.1']
    motifs_cluster_178 = []
    for ma_id in cluster_178:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_178.append(motif)
    motifs_fetched.append(motifs_cluster_178)
    total_count=total_count+len(motifs_cluster_178)

    cluster_179 = ['MA1981.2']
    motifs_cluster_179 = []
    for ma_id in cluster_179:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_179.append(motif)
    motifs_fetched.append(motifs_cluster_179)
    total_count=total_count+len(motifs_cluster_179)

    cluster_180 = ['MA1982.2']
    motifs_cluster_180 = []
    for ma_id in cluster_180:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_180.append(motif)
    motifs_fetched.append(motifs_cluster_180)
    total_count=total_count+len(motifs_cluster_180)

    cluster_181 = ['MA1983.2']
    motifs_cluster_181 = []
    for ma_id in cluster_181:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_181.append(motif)
    motifs_fetched.append(motifs_cluster_181)
    total_count=total_count+len(motifs_cluster_181)

    cluster_182 = ['MA1984.2']
    motifs_cluster_182 = []
    for ma_id in cluster_182:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_182.append(motif)
    motifs_fetched.append(motifs_cluster_182)
    total_count=total_count+len(motifs_cluster_182)

    cluster_183 = ['MA1986.2']
    motifs_cluster_183 = []
    for ma_id in cluster_183:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_183.append(motif)
    motifs_fetched.append(motifs_cluster_183)
    total_count=total_count+len(motifs_cluster_183)

    cluster_184 = ['MA1987.2']
    motifs_cluster_184 = []
    for ma_id in cluster_184:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_184.append(motif)
    motifs_fetched.append(motifs_cluster_184)
    total_count=total_count+len(motifs_cluster_184)

    cluster_185 = ['MA1998.2']
    motifs_cluster_185 = []
    for ma_id in cluster_185:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_185.append(motif)
    motifs_fetched.append(motifs_cluster_185)
    total_count=total_count+len(motifs_cluster_185)

    cluster_186 = ['MA1999.2']
    motifs_cluster_186 = []
    for ma_id in cluster_186:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_186.append(motif)
    motifs_fetched.append(motifs_cluster_186)
    total_count=total_count+len(motifs_cluster_186)

    cluster_187 = ['MA2100.1']
    motifs_cluster_187 = []
    for ma_id in cluster_187:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_187.append(motif)
    motifs_fetched.append(motifs_cluster_187)
    total_count=total_count+len(motifs_cluster_187)

    cluster_188 = ['MA2101.1']
    motifs_cluster_188 = []
    for ma_id in cluster_188:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_188.append(motif)
    motifs_fetched.append(motifs_cluster_188)
    total_count=total_count+len(motifs_cluster_188)

    cluster_189 = ['MA2120.1']
    motifs_cluster_189 = []
    for ma_id in cluster_189:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_189.append(motif)
    motifs_fetched.append(motifs_cluster_189)
    total_count=total_count+len(motifs_cluster_189)

    cluster_190 = ['MA2121.1']
    motifs_cluster_190 = []
    for ma_id in cluster_190:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_190.append(motif)
    motifs_fetched.append(motifs_cluster_190)
    total_count=total_count+len(motifs_cluster_190)

    cluster_191 = ['MA2123.1']
    motifs_cluster_191 = []
    for ma_id in cluster_191:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_191.append(motif)
    motifs_fetched.append(motifs_cluster_191)
    total_count=total_count+len(motifs_cluster_191)

    cluster_192 = ['MA2124.1']
    motifs_cluster_192 = []
    for ma_id in cluster_192:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_192.append(motif)
    motifs_fetched.append(motifs_cluster_192)
    total_count=total_count+len(motifs_cluster_192)

    cluster_193 = ['MA2125.1']
    motifs_cluster_193 = []
    for ma_id in cluster_193:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_193.append(motif)
    motifs_fetched.append(motifs_cluster_193)
    total_count=total_count+len(motifs_cluster_193)

    cluster_194 = ['MA2331.1']
    motifs_cluster_194 = []
    for ma_id in cluster_194:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_194.append(motif)
    motifs_fetched.append(motifs_cluster_194)
    total_count=total_count+len(motifs_cluster_194)

    cluster_195 = ['MA2333.1']
    motifs_cluster_195 = []
    for ma_id in cluster_195:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_195.append(motif)
    motifs_fetched.append(motifs_cluster_195)
    total_count=total_count+len(motifs_cluster_195)

    cluster_196 = ['MA2335.1']
    motifs_cluster_196 = []
    for ma_id in cluster_196:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_196.append(motif)
    motifs_fetched.append(motifs_cluster_196)
    total_count=total_count+len(motifs_cluster_196)

    cluster_197 = ['MA2336.1']
    motifs_cluster_197 = []
    for ma_id in cluster_197:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_197.append(motif)
    motifs_fetched.append(motifs_cluster_197)
    total_count=total_count+len(motifs_cluster_197)

    cluster_198 = ['UN0113.1']
    motifs_cluster_198 = []
    for ma_id in cluster_198:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_198.append(motif)
    motifs_fetched.append(motifs_cluster_198)
    total_count=total_count+len(motifs_cluster_198)

    cluster_199 = ['UN0114.1']
    motifs_cluster_199 = []
    for ma_id in cluster_199:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_199.append(motif)
    motifs_fetched.append(motifs_cluster_199)
    total_count=total_count+len(motifs_cluster_199)

    cluster_200 = ['UN0121.2']
    motifs_cluster_200 = []
    for ma_id in cluster_200:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_200.append(motif)
    motifs_fetched.append(motifs_cluster_200)
    total_count=total_count+len(motifs_cluster_200)

    cluster_201 = ['UN0123.2']
    motifs_cluster_201 = []
    for ma_id in cluster_201:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_201.append(motif)
    motifs_fetched.append(motifs_cluster_201)
    total_count=total_count+len(motifs_cluster_201)

    cluster_202 = ['UN0132.2']
    motifs_cluster_202 = []
    for ma_id in cluster_202:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_202.append(motif)
    motifs_fetched.append(motifs_cluster_202)
    total_count=total_count+len(motifs_cluster_202)

    cluster_203 = ['UN0139.1']
    motifs_cluster_203 = []
    for ma_id in cluster_203:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_203.append(motif)
    motifs_fetched.append(motifs_cluster_203)
    total_count=total_count+len(motifs_cluster_203)

    cluster_204 = ['UN0144.2']
    motifs_cluster_204 = []
    for ma_id in cluster_204:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_204.append(motif)
    motifs_fetched.append(motifs_cluster_204)
    total_count=total_count+len(motifs_cluster_204)

    cluster_205 = ['UN0146.2']
    motifs_cluster_205 = []
    for ma_id in cluster_205:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_205.append(motif)
    motifs_fetched.append(motifs_cluster_205)
    total_count=total_count+len(motifs_cluster_205)

    cluster_206 = ['UN0148.2']
    motifs_cluster_206 = []
    for ma_id in cluster_206:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_206.append(motif)
    motifs_fetched.append(motifs_cluster_206)
    total_count=total_count+len(motifs_cluster_206)

    cluster_207 = ['UN0149.2']
    motifs_cluster_207 = []
    for ma_id in cluster_207:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_207.append(motif)
    motifs_fetched.append(motifs_cluster_207)
    total_count=total_count+len(motifs_cluster_207)

    cluster_208 = ['UN0151.2']
    motifs_cluster_208 = []
    for ma_id in cluster_208:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_208.append(motif)
    motifs_fetched.append(motifs_cluster_208)
    total_count=total_count+len(motifs_cluster_208)

    cluster_209 = ['UN0153.2']
    motifs_cluster_209 = []
    for ma_id in cluster_209:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_209.append(motif)
    motifs_fetched.append(motifs_cluster_209)
    total_count=total_count+len(motifs_cluster_209)

    cluster_210 = ['UN0154.2']
    motifs_cluster_210 = []
    for ma_id in cluster_210:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_210.append(motif)
    motifs_fetched.append(motifs_cluster_210)
    total_count=total_count+len(motifs_cluster_210)

    cluster_211 = ['UN0156.2']
    motifs_cluster_211 = []
    for ma_id in cluster_211:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_211.append(motif)
    motifs_fetched.append(motifs_cluster_211)
    total_count=total_count+len(motifs_cluster_211)

    cluster_212 = ['UN0158.1']
    motifs_cluster_212 = []
    for ma_id in cluster_212:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_212.append(motif)
    motifs_fetched.append(motifs_cluster_212)
    total_count=total_count+len(motifs_cluster_212)

    cluster_213 = ['UN0161.2']
    motifs_cluster_213 = []
    for ma_id in cluster_213:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_213.append(motif)
    motifs_fetched.append(motifs_cluster_213)
    total_count=total_count+len(motifs_cluster_213)

    cluster_214 = ['UN0162.2']
    motifs_cluster_214 = []
    for ma_id in cluster_214:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_214.append(motif)
    motifs_fetched.append(motifs_cluster_214)
    total_count=total_count+len(motifs_cluster_214)

    cluster_215 = ['UN0163.2']
    motifs_cluster_215 = []
    for ma_id in cluster_215:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_215.append(motif)
    motifs_fetched.append(motifs_cluster_215)
    total_count=total_count+len(motifs_cluster_215)

    cluster_216 = ['UN0164.2']
    motifs_cluster_216 = []
    for ma_id in cluster_216:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_216.append(motif)
    motifs_fetched.append(motifs_cluster_216)
    total_count=total_count+len(motifs_cluster_216)

    cluster_217 = ['UN0165.2']
    motifs_cluster_217 = []
    for ma_id in cluster_217:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_217.append(motif)
    motifs_fetched.append(motifs_cluster_217)
    total_count=total_count+len(motifs_cluster_217)

    cluster_218 = ['UN0167.2']
    motifs_cluster_218 = []
    for ma_id in cluster_218:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_218.append(motif)
    motifs_fetched.append(motifs_cluster_218)
    total_count=total_count+len(motifs_cluster_218)

    cluster_219 = ['UN0168.1']
    motifs_cluster_219 = []
    for ma_id in cluster_219:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_219.append(motif)
    motifs_fetched.append(motifs_cluster_219)
    total_count=total_count+len(motifs_cluster_219)

    cluster_220 = ['UN0169.1']
    motifs_cluster_220 = []
    for ma_id in cluster_220:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_220.append(motif)
    motifs_fetched.append(motifs_cluster_220)
    total_count=total_count+len(motifs_cluster_220)

    cluster_221 = ['UN0172.2']
    motifs_cluster_221 = []
    for ma_id in cluster_221:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_221.append(motif)
    motifs_fetched.append(motifs_cluster_221)
    total_count=total_count+len(motifs_cluster_221)

    cluster_222 = ['UN0173.2']
    motifs_cluster_222 = []
    for ma_id in cluster_222:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_222.append(motif)
    motifs_fetched.append(motifs_cluster_222)
    total_count=total_count+len(motifs_cluster_222)

    cluster_223 = ['UN0176.1']
    motifs_cluster_223 = []
    for ma_id in cluster_223:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_223.append(motif)
    motifs_fetched.append(motifs_cluster_223)
    total_count=total_count+len(motifs_cluster_223)

    cluster_224 = ['UN0179.2']
    motifs_cluster_224 = []
    for ma_id in cluster_224:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_224.append(motif)
    motifs_fetched.append(motifs_cluster_224)
    total_count=total_count+len(motifs_cluster_224)

    cluster_225 = ['UN0184.2']
    motifs_cluster_225 = []
    for ma_id in cluster_225:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_225.append(motif)
    motifs_fetched.append(motifs_cluster_225)
    total_count=total_count+len(motifs_cluster_225)

    cluster_226 = ['UN0186.2']
    motifs_cluster_226 = []
    for ma_id in cluster_226:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_226.append(motif)
    motifs_fetched.append(motifs_cluster_226)
    total_count=total_count+len(motifs_cluster_226)

    cluster_227 = ['UN0187.2']
    motifs_cluster_227 = []
    for ma_id in cluster_227:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_227.append(motif)
    motifs_fetched.append(motifs_cluster_227)
    total_count=total_count+len(motifs_cluster_227)

    cluster_228 = ['UN0188.1']
    motifs_cluster_228 = []
    for ma_id in cluster_228:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_228.append(motif)
    motifs_fetched.append(motifs_cluster_228)
    total_count=total_count+len(motifs_cluster_228)

    cluster_229 = ['UN0190.2']
    motifs_cluster_229 = []
    for ma_id in cluster_229:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_229.append(motif)
    motifs_fetched.append(motifs_cluster_229)
    total_count=total_count+len(motifs_cluster_229)

    cluster_230 = ['UN0196.2']
    motifs_cluster_230 = []
    for ma_id in cluster_230:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_230.append(motif)
    motifs_fetched.append(motifs_cluster_230)
    total_count=total_count+len(motifs_cluster_230)

    cluster_231 = ['UN0201.1']
    motifs_cluster_231 = []
    for ma_id in cluster_231:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_231.append(motif)
    motifs_fetched.append(motifs_cluster_231)
    total_count=total_count+len(motifs_cluster_231)

    cluster_232 = ['UN0204.1']
    motifs_cluster_232 = []
    for ma_id in cluster_232:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_232.append(motif)
    motifs_fetched.append(motifs_cluster_232)
    total_count=total_count+len(motifs_cluster_232)

    cluster_233 = ['UN0205.1']
    motifs_cluster_233 = []
    for ma_id in cluster_233:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_233.append(motif)
    motifs_fetched.append(motifs_cluster_233)
    total_count=total_count+len(motifs_cluster_233)

    cluster_234 = ['UN0206.2']
    motifs_cluster_234 = []
    for ma_id in cluster_234:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_234.append(motif)
    motifs_fetched.append(motifs_cluster_234)
    total_count=total_count+len(motifs_cluster_234)

    cluster_235 = ['UN0207.2']
    motifs_cluster_235 = []
    for ma_id in cluster_235:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_235.append(motif)
    motifs_fetched.append(motifs_cluster_235)
    total_count=total_count+len(motifs_cluster_235)

    cluster_236 = ['UN0210.2']
    motifs_cluster_236 = []
    for ma_id in cluster_236:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_236.append(motif)
    motifs_fetched.append(motifs_cluster_236)
    total_count=total_count+len(motifs_cluster_236)

    cluster_237 = ['UN0213.2']
    motifs_cluster_237 = []
    for ma_id in cluster_237:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_237.append(motif)
    motifs_fetched.append(motifs_cluster_237)
    total_count=total_count+len(motifs_cluster_237)

    cluster_238 = ['UN0216.2']
    motifs_cluster_238 = []
    for ma_id in cluster_238:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_238.append(motif)
    motifs_fetched.append(motifs_cluster_238)
    total_count=total_count+len(motifs_cluster_238)

    cluster_239 = ['UN0217.2']
    motifs_cluster_239 = []
    for ma_id in cluster_239:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_239.append(motif)
    motifs_fetched.append(motifs_cluster_239)
    total_count=total_count+len(motifs_cluster_239)

    cluster_240 = ['UN0222.2']
    motifs_cluster_240 = []
    for ma_id in cluster_240:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_240.append(motif)
    motifs_fetched.append(motifs_cluster_240)
    total_count=total_count+len(motifs_cluster_240)

    cluster_241 = ['UN0226.2']
    motifs_cluster_241 = []
    for ma_id in cluster_241:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_241.append(motif)
    motifs_fetched.append(motifs_cluster_241)
    total_count=total_count+len(motifs_cluster_241)

    cluster_242 = ['UN0230.2']
    motifs_cluster_242 = []
    for ma_id in cluster_242:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_242.append(motif)
    motifs_fetched.append(motifs_cluster_242)
    total_count=total_count+len(motifs_cluster_242)

    cluster_243 = ['UN0234.1']
    motifs_cluster_243 = []
    for ma_id in cluster_243:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_243.append(motif)
    motifs_fetched.append(motifs_cluster_243)
    total_count=total_count+len(motifs_cluster_243)

    cluster_244 = ['UN0235.1']
    motifs_cluster_244 = []
    for ma_id in cluster_244:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_244.append(motif)
    motifs_fetched.append(motifs_cluster_244)
    total_count=total_count+len(motifs_cluster_244)

    cluster_245 = ['UN0237.2']
    motifs_cluster_245 = []
    for ma_id in cluster_245:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_245.append(motif)
    motifs_fetched.append(motifs_cluster_245)
    total_count=total_count+len(motifs_cluster_245)

    cluster_246 = ['UN0239.1']
    motifs_cluster_246 = []
    for ma_id in cluster_246:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_246.append(motif)
    motifs_fetched.append(motifs_cluster_246)
    total_count=total_count+len(motifs_cluster_246)

    cluster_247 = ['UN0246.2']
    motifs_cluster_247 = []
    for ma_id in cluster_247:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_247.append(motif)
    motifs_fetched.append(motifs_cluster_247)
    total_count=total_count+len(motifs_cluster_247)

    cluster_248 = ['UN0248.2']
    motifs_cluster_248 = []
    for ma_id in cluster_248:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_248.append(motif)
    motifs_fetched.append(motifs_cluster_248)
    total_count=total_count+len(motifs_cluster_248)

    cluster_249 = ['UN0249.2']
    motifs_cluster_249 = []
    for ma_id in cluster_249:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_249.append(motif)
    motifs_fetched.append(motifs_cluster_249)
    total_count=total_count+len(motifs_cluster_249)

    cluster_250 = ['UN0264.2']
    motifs_cluster_250 = []
    for ma_id in cluster_250:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_250.append(motif)
    motifs_fetched.append(motifs_cluster_250)
    total_count=total_count+len(motifs_cluster_250)

    cluster_251 = ['UN0269.2']
    motifs_cluster_251 = []
    for ma_id in cluster_251:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_251.append(motif)
    motifs_fetched.append(motifs_cluster_251)
    total_count=total_count+len(motifs_cluster_251)

    cluster_252 = ['UN0314.2']
    motifs_cluster_252 = []
    for ma_id in cluster_252:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_252.append(motif)
    motifs_fetched.append(motifs_cluster_252)
    total_count=total_count+len(motifs_cluster_252)

    cluster_253 = ['UN0320.2']
    motifs_cluster_253 = []
    for ma_id in cluster_253:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_253.append(motif)
    motifs_fetched.append(motifs_cluster_253)
    total_count=total_count+len(motifs_cluster_253)

    cluster_254 = ['UN0323.2']
    motifs_cluster_254 = []
    for ma_id in cluster_254:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_254.append(motif)
    motifs_fetched.append(motifs_cluster_254)
    total_count=total_count+len(motifs_cluster_254)

    cluster_255 = ['UN0332.2']
    motifs_cluster_255 = []
    for ma_id in cluster_255:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_255.append(motif)
    motifs_fetched.append(motifs_cluster_255)
    total_count=total_count+len(motifs_cluster_255)

    cluster_256 = ['UN0335.2']
    motifs_cluster_256 = []
    for ma_id in cluster_256:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_256.append(motif)
    motifs_fetched.append(motifs_cluster_256)
    total_count=total_count+len(motifs_cluster_256)

    cluster_257 = ['UN0486.1']
    motifs_cluster_257 = []
    for ma_id in cluster_257:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_257.append(motif)
    motifs_fetched.append(motifs_cluster_257)
    total_count=total_count+len(motifs_cluster_257)

    cluster_258 = ['UN0487.1']
    motifs_cluster_258 = []
    for ma_id in cluster_258:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_258.append(motif)
    motifs_fetched.append(motifs_cluster_258)
    total_count=total_count+len(motifs_cluster_258)

    cluster_259 = ['UN0489.2']
    motifs_cluster_259 = []
    for ma_id in cluster_259:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_259.append(motif)
    motifs_fetched.append(motifs_cluster_259)
    total_count=total_count+len(motifs_cluster_259)

    cluster_260 = ['UN0490.2']
    motifs_cluster_260 = []
    for ma_id in cluster_260:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_260.append(motif)
    motifs_fetched.append(motifs_cluster_260)
    total_count=total_count+len(motifs_cluster_260)

    cluster_261 = ['UN0492.2']
    motifs_cluster_261 = []
    for ma_id in cluster_261:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_261.append(motif)
    motifs_fetched.append(motifs_cluster_261)
    total_count=total_count+len(motifs_cluster_261)

    cluster_262 = ['UN0493.2']
    motifs_cluster_262 = []
    for ma_id in cluster_262:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_262.append(motif)
    motifs_fetched.append(motifs_cluster_262)
    total_count=total_count+len(motifs_cluster_262)

    cluster_263 = ['UN0497.1']
    motifs_cluster_263 = []
    for ma_id in cluster_263:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_263.append(motif)
    motifs_fetched.append(motifs_cluster_263)
    total_count=total_count+len(motifs_cluster_263)

    cluster_264 = ['UN0498.2']
    motifs_cluster_264 = []
    for ma_id in cluster_264:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_264.append(motif)
    motifs_fetched.append(motifs_cluster_264)
    total_count=total_count+len(motifs_cluster_264)

    cluster_265 = ['UN0499.2']
    motifs_cluster_265 = []
    for ma_id in cluster_265:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_265.append(motif)
    motifs_fetched.append(motifs_cluster_265)
    total_count=total_count+len(motifs_cluster_265)

    cluster_266 = ['UN0519.2']
    motifs_cluster_266 = []
    for ma_id in cluster_266:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_266.append(motif)
    motifs_fetched.append(motifs_cluster_266)
    total_count=total_count+len(motifs_cluster_266)

    cluster_267 = ['UN0521.2']
    motifs_cluster_267 = []
    for ma_id in cluster_267:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_267.append(motif)
    motifs_fetched.append(motifs_cluster_267)
    total_count=total_count+len(motifs_cluster_267)

    cluster_268 = ['UN0533.1']
    motifs_cluster_268 = []
    for ma_id in cluster_268:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_268.append(motif)
    motifs_fetched.append(motifs_cluster_268)
    total_count=total_count+len(motifs_cluster_268)

    cluster_269 = ['UN0540.2']
    motifs_cluster_269 = []
    for ma_id in cluster_269:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_269.append(motif)
    motifs_fetched.append(motifs_cluster_269)
    total_count=total_count+len(motifs_cluster_269)

    cluster_270 = ['UN0545.2']
    motifs_cluster_270 = []
    for ma_id in cluster_270:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_270.append(motif)
    motifs_fetched.append(motifs_cluster_270)
    total_count=total_count+len(motifs_cluster_270)

    cluster_271 = ['UN0546.2']
    motifs_cluster_271 = []
    for ma_id in cluster_271:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_271.append(motif)
    motifs_fetched.append(motifs_cluster_271)
    total_count=total_count+len(motifs_cluster_271)

    cluster_272 = ['UN0549.2']
    motifs_cluster_272 = []
    for ma_id in cluster_272:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_272.append(motif)
    motifs_fetched.append(motifs_cluster_272)
    total_count=total_count+len(motifs_cluster_272)

    cluster_273 = ['UN0556.2']
    motifs_cluster_273 = []
    for ma_id in cluster_273:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_273.append(motif)
    motifs_fetched.append(motifs_cluster_273)
    total_count=total_count+len(motifs_cluster_273)

    cluster_274 = ['UN0557.2']
    motifs_cluster_274 = []
    for ma_id in cluster_274:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_274.append(motif)
    motifs_fetched.append(motifs_cluster_274)
    total_count=total_count+len(motifs_cluster_274)

    cluster_275 = ['UN0571.1']
    motifs_cluster_275 = []
    for ma_id in cluster_275:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_275.append(motif)
    motifs_fetched.append(motifs_cluster_275)
    total_count=total_count+len(motifs_cluster_275)

    cluster_276 = ['UN0574.1']
    motifs_cluster_276 = []
    for ma_id in cluster_276:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_276.append(motif)
    motifs_fetched.append(motifs_cluster_276)
    total_count=total_count+len(motifs_cluster_276)

    cluster_277 = ['UN0575.1']
    motifs_cluster_277 = []
    for ma_id in cluster_277:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_277.append(motif)
    motifs_fetched.append(motifs_cluster_277)
    total_count=total_count+len(motifs_cluster_277)

    cluster_278 = ['UN0576.1']
    motifs_cluster_278 = []
    for ma_id in cluster_278:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_278.append(motif)
    motifs_fetched.append(motifs_cluster_278)
    total_count=total_count+len(motifs_cluster_278)

    cluster_279 = ['UN0589.1']
    motifs_cluster_279 = []
    for ma_id in cluster_279:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_279.append(motif)
    motifs_fetched.append(motifs_cluster_279)
    total_count=total_count+len(motifs_cluster_279)

    cluster_280 = ['UN0590.1']
    motifs_cluster_280 = []
    for ma_id in cluster_280:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_280.append(motif)
    motifs_fetched.append(motifs_cluster_280)
    total_count=total_count+len(motifs_cluster_280)

    cluster_281 = ['UN0591.2']
    motifs_cluster_281 = []
    for ma_id in cluster_281:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_281.append(motif)
    motifs_fetched.append(motifs_cluster_281)
    total_count=total_count+len(motifs_cluster_281)

    cluster_282 = ['UN0594.1']
    motifs_cluster_282 = []
    for ma_id in cluster_282:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_282.append(motif)
    motifs_fetched.append(motifs_cluster_282)
    total_count=total_count+len(motifs_cluster_282)

    cluster_283 = ['UN0599.2']
    motifs_cluster_283 = []
    for ma_id in cluster_283:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_283.append(motif)
    motifs_fetched.append(motifs_cluster_283)
    total_count=total_count+len(motifs_cluster_283)

    cluster_284 = ['UN0603.2']
    motifs_cluster_284 = []
    for ma_id in cluster_284:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_284.append(motif)
    motifs_fetched.append(motifs_cluster_284)
    total_count=total_count+len(motifs_cluster_284)

    cluster_285 = ['UN0604.1']
    motifs_cluster_285 = []
    for ma_id in cluster_285:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_285.append(motif)
    motifs_fetched.append(motifs_cluster_285)
    total_count=total_count+len(motifs_cluster_285)

    cluster_286 = ['UN0607.1']
    motifs_cluster_286 = []
    for ma_id in cluster_286:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_286.append(motif)
    motifs_fetched.append(motifs_cluster_286)
    total_count=total_count+len(motifs_cluster_286)

    cluster_287 = ['UN0609.1']
    motifs_cluster_287 = []
    for ma_id in cluster_287:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_287.append(motif)
    motifs_fetched.append(motifs_cluster_287)
    total_count=total_count+len(motifs_cluster_287)

    cluster_288 = ['UN0613.2']
    motifs_cluster_288 = []
    for ma_id in cluster_288:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_288.append(motif)
    motifs_fetched.append(motifs_cluster_288)
    total_count=total_count+len(motifs_cluster_288)

    cluster_289 = ['UN0614.1']
    motifs_cluster_289 = []
    for ma_id in cluster_289:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_289.append(motif)
    motifs_fetched.append(motifs_cluster_289)
    total_count=total_count+len(motifs_cluster_289)

    cluster_290 = ['UN0615.1']
    motifs_cluster_290 = []
    for ma_id in cluster_290:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_290.append(motif)
    motifs_fetched.append(motifs_cluster_290)
    total_count=total_count+len(motifs_cluster_290)

    cluster_291 = ['UN0618.2']
    motifs_cluster_291 = []
    for ma_id in cluster_291:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_291.append(motif)
    motifs_fetched.append(motifs_cluster_291)
    total_count=total_count+len(motifs_cluster_291)

    cluster_292 = ['UN0619.1']
    motifs_cluster_292 = []
    for ma_id in cluster_292:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_292.append(motif)
    motifs_fetched.append(motifs_cluster_292)
    total_count=total_count+len(motifs_cluster_292)

    cluster_293 = ['UN0620.2']
    motifs_cluster_293 = []
    for ma_id in cluster_293:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_293.append(motif)
    motifs_fetched.append(motifs_cluster_293)
    total_count=total_count+len(motifs_cluster_293)

    cluster_294 = ['UN0623.2']
    motifs_cluster_294 = []
    for ma_id in cluster_294:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_294.append(motif)
    motifs_fetched.append(motifs_cluster_294)
    total_count=total_count+len(motifs_cluster_294)

    cluster_295 = ['UN0625.2']
    motifs_cluster_295 = []
    for ma_id in cluster_295:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_295.append(motif)
    motifs_fetched.append(motifs_cluster_295)
    total_count=total_count+len(motifs_cluster_295)

    cluster_296 = ['UN0626.1']
    motifs_cluster_296 = []
    for ma_id in cluster_296:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_296.append(motif)
    motifs_fetched.append(motifs_cluster_296)
    total_count=total_count+len(motifs_cluster_296)

    cluster_297 = ['UN0632.2']
    motifs_cluster_297 = []
    for ma_id in cluster_297:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_297.append(motif)
    motifs_fetched.append(motifs_cluster_297)
    total_count=total_count+len(motifs_cluster_297)

    cluster_298 = ['UN0633.1']
    motifs_cluster_298 = []
    for ma_id in cluster_298:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_298.append(motif)
    motifs_fetched.append(motifs_cluster_298)
    total_count=total_count+len(motifs_cluster_298)

    cluster_299 = ['UN0635.2']
    motifs_cluster_299 = []
    for ma_id in cluster_299:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_299.append(motif)
    motifs_fetched.append(motifs_cluster_299)
    total_count=total_count+len(motifs_cluster_299)

    cluster_300 = ['UN0636.2']
    motifs_cluster_300 = []
    for ma_id in cluster_300:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_300.append(motif)
    motifs_fetched.append(motifs_cluster_300)
    total_count=total_count+len(motifs_cluster_300)

    cluster_301 = ['UN0637.2']
    motifs_cluster_301 = []
    for ma_id in cluster_301:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_301.append(motif)
    motifs_fetched.append(motifs_cluster_301)
    total_count=total_count+len(motifs_cluster_301)

    cluster_302 = ['UN0640.2']
    motifs_cluster_302 = []
    for ma_id in cluster_302:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_302.append(motif)
    motifs_fetched.append(motifs_cluster_302)
    total_count=total_count+len(motifs_cluster_302)

    cluster_303 = ['UN0641.2']
    motifs_cluster_303 = []
    for ma_id in cluster_303:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_303.append(motif)
    motifs_fetched.append(motifs_cluster_303)
    total_count=total_count+len(motifs_cluster_303)

    cluster_304 = ['UN0643.1']
    motifs_cluster_304 = []
    for ma_id in cluster_304:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_304.append(motif)
    motifs_fetched.append(motifs_cluster_304)
    total_count=total_count+len(motifs_cluster_304)

    cluster_305 = ['UN0644.2']
    motifs_cluster_305 = []
    for ma_id in cluster_305:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_305.append(motif)
    motifs_fetched.append(motifs_cluster_305)
    total_count=total_count+len(motifs_cluster_305)

    cluster_306 = ['UN0646.2']
    motifs_cluster_306 = []
    for ma_id in cluster_306:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_306.append(motif)
    motifs_fetched.append(motifs_cluster_306)
    total_count=total_count+len(motifs_cluster_306)

    cluster_307 = ['UN0652.2']
    motifs_cluster_307 = []
    for ma_id in cluster_307:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_307.append(motif)
    motifs_fetched.append(motifs_cluster_307)
    total_count=total_count+len(motifs_cluster_307)

    cluster_308 = ['UN0656.2']
    motifs_cluster_308 = []
    for ma_id in cluster_308:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_308.append(motif)
    motifs_fetched.append(motifs_cluster_308)
    total_count=total_count+len(motifs_cluster_308)

    cluster_309 = ['UN0657.2']
    motifs_cluster_309 = []
    for ma_id in cluster_309:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_309.append(motif)
    motifs_fetched.append(motifs_cluster_309)
    total_count=total_count+len(motifs_cluster_309)

    cluster_310 = ['UN0660.2']
    motifs_cluster_310 = []
    for ma_id in cluster_310:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_310.append(motif)
    motifs_fetched.append(motifs_cluster_310)
    total_count=total_count+len(motifs_cluster_310)

    cluster_311 = ['UN0661.2']
    motifs_cluster_311 = []
    for ma_id in cluster_311:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_311.append(motif)
    motifs_fetched.append(motifs_cluster_311)
    total_count=total_count+len(motifs_cluster_311)

    cluster_312 = ['UN0664.2']
    motifs_cluster_312 = []
    for ma_id in cluster_312:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_312.append(motif)
    motifs_fetched.append(motifs_cluster_312)
    total_count=total_count+len(motifs_cluster_312)

    cluster_313 = ['UN0678.2']
    motifs_cluster_313 = []
    for ma_id in cluster_313:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_313.append(motif)
    motifs_fetched.append(motifs_cluster_313)
    total_count=total_count+len(motifs_cluster_313)

    cluster_314 = ['UN0679.1']
    motifs_cluster_314 = []
    for ma_id in cluster_314:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_314.append(motif)
    motifs_fetched.append(motifs_cluster_314)
    total_count=total_count+len(motifs_cluster_314)

    cluster_315 = ['UN0805.1']
    motifs_cluster_315 = []
    for ma_id in cluster_315:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_315.append(motif)
    motifs_fetched.append(motifs_cluster_315)
    total_count=total_count+len(motifs_cluster_315)

    cluster_316 = ['UN0810.1']
    motifs_cluster_316 = []
    for ma_id in cluster_316:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_316.append(motif)
    motifs_fetched.append(motifs_cluster_316)
    total_count=total_count+len(motifs_cluster_316)

    cluster_317 = ['UN0813.1']
    motifs_cluster_317 = []
    for ma_id in cluster_317:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_317.append(motif)
    motifs_fetched.append(motifs_cluster_317)
    total_count=total_count+len(motifs_cluster_317)

    cluster_318 = ['UN0815.1']
    motifs_cluster_318 = []
    for ma_id in cluster_318:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_318.append(motif)
    motifs_fetched.append(motifs_cluster_318)
    total_count=total_count+len(motifs_cluster_318)

    cluster_319 = ['UN0816.1']
    motifs_cluster_319 = []
    for ma_id in cluster_319:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_319.append(motif)
    motifs_fetched.append(motifs_cluster_319)
    total_count=total_count+len(motifs_cluster_319)

    cluster_320 = ['UN0817.1']
    motifs_cluster_320 = []
    for ma_id in cluster_320:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_320.append(motif)
    motifs_fetched.append(motifs_cluster_320)
    total_count=total_count+len(motifs_cluster_320)

    cluster_321 = ['UN0818.1']
    motifs_cluster_321 = []
    for ma_id in cluster_321:
        motif = jdb.fetch_motif_by_id(ma_id)
        motifs_cluster_321.append(motif)
    motifs_fetched.append(motifs_cluster_321)
    total_count=total_count+len(motifs_cluster_321)

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
    
    self.total_cluster_count=321
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


# db = DNABindingMotifs()
# sth1 = db.mmm[list(db.mmm.keys())[0]]
# sth2 = db.mmm[list(db.mmm.keys())[0]]
# print(sth1.pssm)
# sth1.pseudocounts = motifs.jaspar.calculate_pseudocounts(sth1)
# sth2.pseudocounts = {'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
# print(sth1.pssm)
# print(sth1.pssm.dist_pearson(sth2.pssm))
