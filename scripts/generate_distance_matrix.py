"""."""
from __future__ import print_function
import pandas as pd
from string_kernel import io
from string_kernel.core.src import string_kernel, sum_string_kernel

df = io.pdb_to_df('/home/fede/projects_local/fais/data_CLL_2015/models')

# df_tmp = pd.read_csv("/home/fede/Dropbox/projects/Franco_Fabio_Marcat/tab_final_merged_newid_mutation.csv")
df_tmp = pd.read_csv("/home/fede/Dropbox/projects/Franco_Fabio_Marcat/tab_01-09-16-nodup.csv")

# new_ids = []
# folders = map(lambda x: x, list(df.index))
# for i in folders:
#     row = list(df_tmp[df_tmp['Folder PAOLO'] == i]['new_id_tm']) or \
#         list(df_tmp[df_tmp['ID TM matrice'] == i[:-4]]['new_id_tm']) or \
#         list(df_tmp[df_tmp['key'] == i[:-4]]['new_id_tm']) or \
#         (print("no correspondence for ", i) or [i[:-4] + '___'])
#     new_ids.append(row[0])
#     if row[0][:4] == 'G034':
#         print(i, row[0])

# sort df according to df_tmp['key']
new_ids = []

loops = ['LCDR1', 'LCDR2', 'LCDR3', 'HCDR1', 'HCDR2', 'HCDR3']
new_df = pd.DataFrame(columns=loops)
for i, fold in zip(df_tmp['key'], df_tmp['Folder PAOLO']):
    new_df = new_df.append(df.loc[fold], ignore_index=True)
    new_ids.append(i)

min_kn = 1
max_kn = 5
lamda = .75
normalize = 1
for i in new_df.columns:
    sum_string_kernel.sum_string_kernel(
        list(new_df[i]),
        filename='{}_kn{}-{}_l{}_norm{}_nodup.csv'.format(i, min_kn, max_kn, lamda, normalize),
        min_kn=min_kn, max_kn=max_kn, lamda=lamda, labels=new_ids,
        normalize=normalize
    )
