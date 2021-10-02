from qiime2.plugins import picrust2
from qiime2 import Artifact

seq_tab = Artifact.import_data("FeatureData[Sequence]", "data/gevers_seq.fa")
otu_tab = Artifact.import_data("FeatureTable[Frequency]", "data/gevers_biom.biom", view_type= 'BIOMV100Format')

picrust2.methods.full_pipeline(table = otu_tab, seq = seq_tab, hsp_method = "pic", threads = 2, placement_tool = "sepp")
