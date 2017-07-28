"""Utilities for Franco IG analysis."""
import pandas as pd
import matplotlib.pyplot as plt

from collections import namedtuple
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

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


def dfcluster_to_namedtuple(cluster, df, verbose=False, col_os="OS",
                            col_dead="Alive/Dead"):
    NamedTuple = namedtuple(
        'NamedTuple', ['id', 'os', 'dead', 'row'], verbose=False)
    myclust_filled = []
    for _id in cluster:
        # row = df[df['new_id_tm'] == _id]
        row = df[df['key'] == _id]
        try:
            os = list(row[col_os])[0]
            float(os)
            dead = list(row[col_dead])[0]
            if os != '':
                myclust_filled.append(NamedTuple(_id, os, dead, row))
        except:
            if verbose:
                print("No correspondence for", _id)
    return sorted(myclust_filled, key=lambda x: x.os)

def print_kaplans_from_clusters(cluster_list, data, ax=None, df=None,
                                alpha=1e-14, col_os="OS",
                                col_dead="Alive/Dead",
                                label=None, verbose=False, **kwargs):
    if cluster_list == []:
        return
    if not isinstance(cluster_list, list):
        cluster_list = [cluster_list]
    kmf = KaplanMeierFitter(alpha=alpha)
    __df = pd.DataFrame(dfcluster_to_namedtuple(
        cluster_list[0], data, verbose=verbose,
        col_os=col_os,
        col_dead=col_dead)) if df is None else df
    if len(__df) == 0:
        return ax

    try:
        kmf.fit(__df['os'], event_observed=__df['dead'], label=label)
        label += '_' + str(kmf.median_)
        kmf.fit(__df['os'], event_observed=__df['dead'], label=label)
        if ax is not None:
            kwargs['ax'] = ax
        ax = kmf.plot(**kwargs)

        for i, clust in enumerate(cluster_list[1:]):
            assert False
            __df = pd.DataFrame(dfcluster_to_namedtuple(
                clust, data, col_os=col_os,
                col_dead=col_dead))
            kmf.fit(__df['os'], event_observed=__df['dead'])
            ax = kmf.plot(ax=ax)
    except Exception as e:
        print("Error for", __df)
    return ax, __df


def kaplan(dataframe, col_condition, col_os='os', col_dead='dead', ax=None,
           alpha=1e-14):
    """Plot a kaplan meier."""
    values = sorted(dataframe[col_condition].unique())
    if ax is None:
        _, ax = plt.subplots()

    kmf = KaplanMeierFitter(alpha=alpha)
    cond = dataframe[col_os].notnull() & dataframe[col_dead].notnull()
    for label in values:
        cond_i = dataframe[col_condition] == label
        df_i = dataframe[cond & cond_i]
        T = df_i[col_os].astype(float)
        E = df_i[col_dead].astype(float)
        kmf.fit(T, E, label=label)
        label += '_' + str(kmf.median_)
        kmf.fit(T, E, label=label)  # this is for labeling
        kmf.plot(ax=ax, show_censors=True)

    ax.set_ylabel("survival probability")
    if len(values) == 2:
        results = logrank_test(
            dataframe[cond & (dataframe[col_condition] == values[0])][col_os],
            dataframe[cond & (dataframe[col_condition] == values[1])][col_os],
            dataframe[cond & (dataframe[col_condition] == values[0])][col_dead],
            dataframe[cond & (dataframe[col_condition] == values[1])][col_dead],
            alpha=.99)
        ax.text(0.1, -0.05, "logrank_test: %.2e" % results.p_value,
                fontsize='medium')
    return ax
