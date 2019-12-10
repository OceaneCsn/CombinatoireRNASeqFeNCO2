import sys
import os
import numpy as np
import pandas as pd


#----------------------------------- Merge many htseq count output in one matrix wit headers -----------------------------------#


#Nom des fichiers sra a traiter
samples = os.listdir(sys.argv[1])
print("Going to merge : \n ",samples)

df = pd.read_csv(sys.argv[1]+"/"+samples[0], sep = '\t', index_col = 0, names =  ["Gene", samples[0].split('.')[0]])
print(df)
for sample in samples[1:]:
    print(sample)
    tmp = pd.read_csv(sys.argv[1]+"/"+sample, sep = '\t', index_col = 0, names = ["Gene", sample.split('.')[0]])    
    df = df.merge(tmp, right_index = True, left_index = True)

print(df)
df.to_csv(sys.argv[2], index=True, header=True)