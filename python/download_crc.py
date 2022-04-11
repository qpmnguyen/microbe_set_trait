from tqdm.autonotebook import tqdm
import pandas as pd 
import os
import matplotlib.pyplot as plt
import subprocess
from multiprocessing import Pool
import time



def download_sample(idx):
    failed = []
    for j in range(0,2):
        url = manifest.iloc[idx].fastq_ftp.split(";")[j]
        sname = sname = manifest.iloc[idx].sample_title.split(" ")[-1]
        print(j)
        if j == 0:
            fname = sname + "_R1_001.fastq.gz"
        elif j == 1:
            fname = sname + "_R2_001.fastq.gz"
        print(fname)
        fullname = crc_path + "/" + fname
        cmd = "wget {} -O {}".format(url, crc_path + "/" + fname)
        if os.path.exists(fullname):
            pass
        else:
            print("Downloading {}".format(fname))
            if idx % 10 == 0 & idx >= 10:
                time.sleep(10)
            subprocess.run(args=["wget", url, "-O", fullname, "--quiet"], stdout=subprocess.DEVNULL)
            if os.path.getsize(fullname) == 0:
                print("For some reason this is not downloading, retrying...")
                subprocess.run(args=["wget", url, "-O", fullname, "--quiet"], stdout=subprocess.DEVNULL)
                if os.path.getsize(fullname) == 0:
                    print("This file is dud")
                    failed.append(url)
                    pass
    print(failed)
    return(failed)



if __name__ == "__main__":
    basepath = "/dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles/data"
    asp_key = "/dartfs-hpc/rc/home/k/f00345k/.aspera/connect/etc/asperaweb_id_dsa.openssh"
    asp_cmd = "ascp -k 1 -QT -l 300m -P33001 -i" 
    manifest = pd.read_csv("crc_16s.tsv", sep="\t")
    manifest = manifest[manifest.library_strategy == "AMPLICON"]
    manifest = manifest[manifest.groupby("sample_title")['read_count'].transform('max') == manifest['read_count']]
    manifest = manifest.reset_index().drop('index', axis = 1)
    crc_path = basepath + "/crc_16s"
    for i in range(0, manifest.shape[0]):
        download_sample(i)
    
    #pool = Pool(5)
    #pool.map(download_sample, range(0, manifest.shape[0]))
    #pool.close()
    #pool.join()