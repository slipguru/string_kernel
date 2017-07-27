"""Utilities for Franco IG analysis."""
import pandas as pd


# old functions, needed to calculate kaplan-meier manually
def namedtuple_to_data_kaplan(cluster, verbose=False):
    if cluster == []:
        return None
    si = 1.
    KM_df = pd.DataFrame(columns=['i', 'ti', 'di', 'ni', 'si'])
    ni = len(cluster)
    cumulative_losses = 0

    # create dict; this allows to group data with same OS
    ti_dict = {}
    for x in cluster:
        ti_dict.setdefault(x.os, []).append(x)

    # order data for OS, then build dataframe for kaplan-meier
    for i, (k, v) in enumerate(sorted(ti_dict.iteritems(), key=lambda t: t[0])):
        deaths = len([_ for _ in v if _.dead == 1.])
        losts = len([_ for _ in v if _.dead == 0.])
        if deaths == 0 and losts > 0:
            # this is the case in which the only information is on information loss;
            # not create this row, but remember the patient loss for the next row
            cumulative_losses = losts
            continue
        ni -= cumulative_losses  # at previous time someone was lost, remove them
        cumulative_losses = 0
        KM_df = KM_df.append({'i': i,'ti': k,'di': deaths,'ni':ni, 'si':si}, ignore_index=True)
        si *= (ni - deaths) / ni
        ni -= deaths - losts
    if verbose: print(KM_df)
    return KM_df

def print_kaplans(df_list):
    if isinstance(df_list, pd.DataFrame): # with a single dataframe
        df_list = [df_list]
    if df_list == []:
        print("No data")

    for df in df_list:
        if df is None or len(df) == 0: continue
        plt.step(df['ti'], df['si'])
    plt.ylim((0,1.1));

# print_kaplans([namedtuple_to_data_kaplan(dfcluster_to_namedtuple(x)) for x in all_clusts])
