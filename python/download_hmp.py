import pandas as pd 
import os 
import subprocess
from multiprocessing import Pool

def download_hmp(idx):
    url = url_file.urls[idx]
    print(url)
    sample_name = url_file.samples[idx] + ".tar.bz2"
    fullname = bpath + "/" + sample_name
    print("Downloading {}".format(fullname))
    if os.path.exists(fullname):
        pass
    else: 
        subprocess.run(args=["wget", url, "-O", fullname, "--quiet"], 
                       stdout = subprocess.DEVNULL)
    

    
if __name__ == "__main__":
    bpath = "/dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles/data/hmp_16s"
    url_file = pd.read_csv("hmp_urls.csv", index_col=0)
    print(url_file.head())
    pool = Pool(5)
    pool.map(download_hmp, range(0, url_file.shape[0]))
    pool.close()
    pool.join()
    