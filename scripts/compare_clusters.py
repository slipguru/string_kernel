from __future__ import print_function
import sys
import numpy as np
import pandas as pd

from icing.plotting import silhouette
from icing.plotting.silhouette import best_intersection

def load_cluster_from_compressed(filename, col, name_id_col='id'):
    df = pd.read_csv(filename, usecols=(name_id_col, col))

    ids = np.array(df[name_id_col])
    cluster_labels = np.array(df[col])
    clusters = {}
    for i in np.unique(cluster_labels):
        clusters[i] = ids[np.where(cluster_labels == i)]

    return clusters


def overlapping(clusts1, clusts2):
    stability_num = 0.  # numerator
    stability_den = 0.  # denominator
    nclusts = len(clusts1)
    orig_clusts_2 = clusts1.copy()
    orig_clusts = []
    for _ in clusts2.itervalues():
        res, clusts1, best_set = best_intersection(_, clusts1)
        n_unknown = len([xx for xx in best_set if xx.endswith("___")])
        print("{1:.2f}"  # (K {2:.2f}, L {3:.2f})"
              .format(_, res[0],
                      (len(best_set) in (0, n_unknown) and -1) or
                      (len([xx for xx in best_set if xx.endswith("_K") or
                      xx.endswith("_")]) - n_unknown) * 100. /
                      (len(best_set) - n_unknown),
                      (len(best_set) in (0, n_unknown) and -1) or
                      (len([xx for xx in best_set if xx.endswith("_L") or
                      xx.endswith("_")]) - n_unknown) * 100. /
                      (len(best_set) - n_unknown)
                      ), end=' ')
        stability_num += res[1]
        stability_den += res[2]
        # print(res, end=' ')
        orig_clusts.append(best_set)

    km, lm = [], []
    for c in orig_clusts_2.itervalues():
        if len(c) > 0:
            n_unknown = len([xx for xx in c if xx.endswith("_")])
            if n_unknown < len(c):
                k_res = (len([xx for xx in c if xx.endswith("_K") or
                         xx.endswith("_")])-n_unknown) * 100. / (len(c)-n_unknown)

                l_res = (len([xx for xx in c if xx.endswith("_L") or
                         xx.endswith("_")])-n_unknown) * 100. / (len(c)-n_unknown)
                if k_res > l_res:
                    km.append(k_res)
                elif k_res < l_res:
                    lm.append(l_res)
                else:
                    km.append(k_res)
                    lm.append(l_res)

    mm, um = [], []
    for c in orig_clusts_2.itervalues():
        if len(c) > 0:
            n_unknown = len([xx for xx in c if xx.endswith('__')])
            if n_unknown < len(c):
                m_res = (len([xx for xx in c if '_M_' in xx or
                         xx.endswith('__')])-n_unknown) * 100. / (len(c)-n_unknown)

                u_res = (len([xx for xx in c if '_U_' in xx or
                         xx.endswith('__')])-n_unknown) * 100. / (len(c)-n_unknown)
                if m_res > u_res:
                    mm.append(m_res)
                elif m_res < u_res:
                    um.append(u_res)
                else:
                    mm.append(m_res)
                    um.append(u_res)

    stability = stability_num / stability_den if stability_den != 0 else 0
    print("\nstability: {:.3f}, {} clusts (light), {} nclusts -- "
          "[{:.2f}%] k[{:.2f}] l[{:.2f}] m[{:.2f}] u[{:.2f}]"
          .format(stability, len(clusts2), nclusts,
                  stability * 100., np.mean(km), np.mean(lm), np.mean(mm), np.mean(um)))


def main():
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    col1 = sys.argv[3]
    col2 = sys.argv[4]

    clust1 = load_cluster_from_compressed(filename1, col1)
    clust2 = load_cluster_from_compressed(filename2, col2)
    overlapping(clust1, clust2)


if __name__ == '__main__':
    main()
