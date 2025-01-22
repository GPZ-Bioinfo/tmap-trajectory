import markov_clustering as mc
import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import MinMaxScaler


def upgma(u: np.array, v: np.array, metric='euclidean') -> float:
    """
    UPGMA algorithm to calculate the distance between cluster u and cluster v.
    For all points and where and are the cardinalities of clusters and , respectively.
    
    d(u, v) = \Sigma_i,j(d(u[i], v[j]) / |u|*|v|)
    
    Refers to average linkage https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage

    Args:
        u (np.array): _description_
        v (np.array): _description_

    Returns:
        float: _description_
    """
    if u.ndim == 1:
        # transform row vector to column vector, aka. reshape to (-1, 1)
        u = u[:, None]
    
    if v.ndim == 1:
        v = v[:, None]

    assert u.shape[1] == v.shape[1]

    x = np.concatenate((u, v))
    mat = squareform(pdist(x, metric=metric))
    dis = (mat[:u.shape[0], -v.shape[0]:] / (u.shape[0]*v.shape[0])).sum()
    return dis


def mcl_cluster(g: nx.Graph, expansion: int = None, inflation: float = None) -> dict:
    mat = csr_matrix(nx.to_scipy_sparse_array(g))

    def clustering(**kwargs):
        res = mc.run_mcl(mat, **kwargs)
        clus = mc.get_clusters(res)
        qual = mc.modularity(matrix=res, clusters=clus)
        print(f'{kwargs}\tmodularity: {qual}')
        return clus, qual

    if expansion is None:
        exps = np.arange(2, 7)
        qs = list(map(lambda x: clustering(expansion=x)[1], exps))
        expansion = exps[np.argmax(qs)]
    
    if inflation is None:
        infs = np.arange(1.5, 3, 0.1)
        qs = list(map(lambda x: clustering(expansion=expansion, inflation=x)[1], infs))
        inflation = infs[np.argmax(qs)]

    clus, qual = clustering(expansion=expansion, inflation=inflation)
    return {k: f'cluster {v}' for v, l in enumerate(clus) for k in l}, {'expansion': expansion, 'inflation': inflation, 'modularity': qual}


if __name__ == '__main__':
    import itertools as its
    import pandas as pd
    import pickle
    import plotly.graph_objs as go
    import plotly.express as px
    from collections import Counter
    from plotly.subplots import make_subplots
    from .simplify_digraph import export_graph, vis_digraph

    age = pd.read_csv('data/metadata.csv', index_col=0)['age_months']
    with open('data/tmap_graph_r35_eps99_overlap70.pkl', 'rb') as f:
        graph = pickle.load(f)
    
    dis = np.array([
        upgma(u=age.loc[graph.node2sample(src)].values, 
              v=age.loc[graph.node2sample(dst)].values) \
        for src, dst in graph.edges])
    mms = MinMaxScaler()
    dis_scaled = mms.fit_transform(dis[:, None]).flatten()

    node_size = MinMaxScaler((10, 40)).fit_transform(np.array([graph.nodes[_].get('size') for _ in list(graph)]).reshape(-1, 1)).flatten()
    pos = dict(zip(graph.nodes, graph.nodePos))

    # 1. unweighted MCL
    nodelink = nx.node_link_data(graph)
    g = nx.node_link_graph(nodelink)
    clus_uw, params = mcl_cluster(g)
    nx.set_node_attributes(
        g,
        values=clus_uw,
        name='MCL clusters (unweighted)'
    )
    export_graph(g, 'data/mcl_graph_unweighted')
    node_val = [clus_uw.get(_) for _ in list(g)]
    fig = vis_digraph(g, node_val, node_size, pos, arrows=False)
    fig.layout.title = f"MCL (unweighted)<BR>expansion: {params.get('expansion'):.2f}, inflation: {params.get('inflation'):.2f}, modularity: {params.get('modularity'):.2f}"
    fig.write_html('data/fig_tmap_graph_colored_mcl_unweighted.html')

    # 2. age scaled MCL
    g = nx.node_link_graph(nodelink)
    nx.set_edge_attributes(
        g,
        values={link: {'weight': wt} for link, wt in zip(graph.edges, 1-dis_scaled)}
    )
    clus_wt, params = mcl_cluster(g)
    nx.set_node_attributes(
        g,
        values=clus_wt,
        name='MCL clusters (age scaled)'
    )
    export_graph(g, 'data/mcl_graph_age_scaled')
    node_val = [clus_wt.get(_) for _ in list(g)]
    fig = vis_digraph(g, node_val, node_size, pos, arrows=False)
    fig.layout.title = f"MCL (age scaled)<BR>expansion: {params.get('expansion'):.2f}, inflation: {params.get('inflation'):.2f}, modularity: {params.get('modularity'):.2f}"
    fig.write_html('data/fig_tmap_graph_colored_mcl_age_scaled.html')

    # 3. Age inverse MCL
    g = nx.node_link_graph(nodelink)
    nx.set_edge_attributes(
        g,
        values={link: {'weight': wt} for link, wt in zip(graph.edges, 1.0/(1 + dis))}
    )
    clus_wt_inv, params = mcl_cluster(g)
    nx.set_node_attributes(
        g,
        values=clus_wt_inv,
        name='MCL clusters (age inverse)'
    )
    export_graph(g, 'data/mcl_graph_age_inverse')
    node_val = [clus_wt_inv.get(_) for _ in list(g)]
    fig = vis_digraph(g, node_val, node_size, pos, arrows=False)
    fig.layout.title = f"MCL (age inverse)<BR>expansion: {params.get('expansion'):.2f}, inflation: {params.get('inflation'):.2f}, modularity: {params.get('modularity'):.2f}"
    fig.write_html('data/fig_tmap_graph_colored_mcl_age_inverse.html')

    # next, measure time labels on clusters
    df = pd.DataFrame.from_dict(
        {'Age': {_: age.loc[graph.node2sample(_)].mean() for _ in list(g)}, 'MCL clusters (unweighted)':  clus_uw, 'MCL clusters (age scaled)': clus_wt, 'MCL clusters (age inverse)': clus_wt_inv},
        orient='columns')
    df.to_csv('data/mcl_time_label_consistency.csv', index_label='NodeID')
    figs = make_subplots(rows=3, cols=1, subplot_titles=['MCL clusters (unweighted)', 'MCL clusters (age scaled)', 'MCL clusters (age inverse)'])
    
    for i, nm in enumerate(['MCL clusters (unweighted)', 'MCL clusters (age scaled)', 'MCL clusters (age inverse)']):
        cmap = {cat: col for (cat, cnt), col in zip(Counter(df[nm]).most_common(), its.cycle(px.colors.qualitative.Plotly[1:]+[px.colors.qualitative.Plotly[0]]))}
        clus_ord = df.groupby(nm).Age.mean().sort_values().index[::-1]
        clus_desc = df.groupby(nm).Age.describe().loc[clus_ord]
        [figs.add_trace(go.Violin(
                x=df[nm][df[nm] == cat],
                y=df['Age'][df[nm] == cat],
                name=cat,
                line_color=cmap.get(cat),
                box_visible=True,
                meanline_visible=True,
            ), row=i+1, col=1) \
            for cat in clus_ord]

    figs.layout.title = "Time label consistency across MCL clusters"
    figs.layout.width = 1200
    figs.layout.height = 1000
    figs.layout.showlegend = False
    figs.write_html('data/fig_mcl_time_label_consistency.html')
