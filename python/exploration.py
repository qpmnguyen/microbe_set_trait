from Bio import SeqIO 
import gzip 

with gzip.open("/Users/quangnguyen/Downloads/ena_files/ERR1368879/1939.100001.fastq.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        print(record.id)