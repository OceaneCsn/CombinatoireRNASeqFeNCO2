import sys
import os
import pandas as pd
import numpy as np

files = os.listdir("./Logs")
df = pd.DataFrame(np.zeros([len(files), 4], dtype=float), columns = ["ID","UM", "MM", "MappedReads"])

for i,f in enumerate(files):
    MappedReads = 0
    print('name : ', f.split('R')[1].split('Log')[0])
    ID = int(f.split('R')[1].split('Log')[0])
    df.loc[i, "ID"] = ID
    log = open("./Logs/"+f, "r")
    for line in log.readlines():
        if line != "":
            if "Number of reads mapped to multiple loci
            
            " in line:
                df.loc[i, "MappedReads"]+= float(line.split('\t')[1].split('\n')[0])
            elif "of reads mapped to multiple loci"in line:
                df.loc[i, "MM"]+= float(line.split('\t')[1].split('%')[0])
            if "Uniquely mapped reads number" in line:
                df.loc[i, "MappedReads"]+= float(line.split('\t')[1].split('\n')[0])
            elif "Uniquely mapped reads"in line:
                df.loc[i, "UM"]+= float(line.split('\t')[1].split('%')[0])

print(df)
df.to_csv("./LogsData.csv", index=False, header=True)