"""."""
import pandas as pd
from string_kernel import io
from string_kernel.core.src import string_kernel, sum_string_kernel

df = io.pdb_to_df('/home/fede/Downloads/data_CLL_2015/data_CLL_2015/models_rf (copy)')

df_tmp = pd.read_csv("/home/fede/Dropbox/projects/Franco_Fabio_Marcat/conversioni_ID/tab_final_merged_newid_mutation.csv", usecols=("Folder PAOLO", "new_id_tm"))

new_ids = []
folders = map(lambda x: x, list(df.index))
for i in folders:
    row = df_tmp[df_tmp['Folder PAOLO'] == i]
    try:
        new_ids.append(list(row['new_id_tm'])[0])
    except:
        print("no correspondece for", i)
        new_ids.append(i[:-4] + '___')

for i in df.columns:
    sum_string_kernel.sum_string_kernel(
        list(df[i]), filename='{}_kn1-3_l0.5.csv'.format(i), min_kn=1, max_kn=3,
        lamda=0.5, labels=new_ids
    )
