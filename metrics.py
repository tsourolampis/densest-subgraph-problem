import math
from networkx.algorithms import bipartite

def fe_subgraph(G, S):
    '''fe(G[S])= e(S)/(|S| choose 2)'''
    H = G.subgraph(S)
    n = len(S)
    return 2 * H.number_of_edges() / float(n * (n - 1))


def fe_bipartite_subgraph(G, S=None, T=None):
    '''fe(G[S \cup T ])= e(S,T)/|S||T|  '''
    if (bipartite.is_bipartite_node_set(G, S) == False):
        return -1  # sentinel value, o/w value is between 0,1

    s, S, t, T = len(S), set(S), len(T), set(T)
    H = G.subgraph(list(S | T))
    return H.number_of_edges() / (s * t)


def rho_subgraph(G, S):
    ''' degree density e(S)/|S| '''
    H = G.subgraph(S)
    return H.number_of_edges() / float(H.number_of_nodes())


def rho_bipartite_subgraph(G, S, T):
    ''' degree density e(S,T)/sqrt(|S||T|) '''
    if (bipartite.is_bipartite_node_set(G, S) == False):
        return -1
    s, S, t, T = len(S), set(S), len(T), set(T)
    H = G.subgraph(list(S | T))
    return H.number_of_edges() / (math.sqrt(s * t))


def edge_surplus_subgraph(G, S, alpha=0.1, gamma=2):
    """
    Returns the edge surplus of the graph edge surplus
    e(S) - alpha*(|S| choose 2) for gamma=2 o/w
    e(S) - alpha* |S|^gamma for gamma=2 o/w
    Suggested values of gamma fall in the range [1,2]
    """
    H = G.subgraph(S)
    n = len(S)
    if gamma == 2:
        return H.number_of_edges() - alpha * (n * (n - 1)) // 2
    else:
        return H.number_of_edges() - alpha * math.pow(n, gamma)


def edge_surplus_bipartite_subgraph(G, S, T, alpha=0.1, gamma=1):
    """
    Returns the edge surplus of the graph edge surplus
    e(S,T) - gamma*|S||T|
    """
    if (bipartite.is_bipartite_node_set(G, S) == False):
        return -1
    s, S, t, T = len(S), set(S), len(T), set(T)
    H = G.subgraph(list(S | T))
    return H.number_of_edges() - alpha * math.pow(s * t, gamma)


def fe(G):
    return fe_subgraph(G, G.nodes())


def fe_bipartite(G):
    if (bipartite.is_bipartite(G) == False):
        return -1
    L, R = bipartite.sets(G)
    return fe_bipartite_subgraph(G, L, R)


def rho(G):
    return rho_subgraph(G, G.nodes())


def rho_bipartite(G):
    if (bipartite.is_bipartite(G) == False):
        return -1
    L, R = bipartite.sets(G)
    return rho_bipartite_subgraph(G, L, R)


def edge_surplus(G, alpha=0.1, gamma=2):
    """
    Returns the edge surplus of the graph edge surplus
    e(S) - gamma*(|S| choose 2)
    """
    return edge_surplus_subgraph(G, G.nodes(), alpha, gamma)


def edge_surplus_bipartite(G, alpha=0.1, gamma=1):
    """
    Returns the edge surplus of a bipartite graph G
    """
    if (bipartite.is_bipartite(G) == False):
        return -1
    L, R = bipartite.sets(G)
    return edge_surplus_bipartite_subgraph(G, L, R, alpha, gamma)


def denstats(G, bipartiteflag=None, alpha=0.1, gamma=2):
    if (bipartiteflag == None):  # avoid testing if we already know, if G is bipartite or not
        bipartiteflag = bipartite.is_bipartite(G)

    if (bipartiteflag == False):
        d, s, fedge = rho(G), edge_surplus(G, alpha, gamma), fe(G)
        print('Non-bipartite graph statistics')
        print('n,m :', G.number_of_nodes(), G.number_of_edges())
        print('rho(G): ', d)
        print('edge-surplus(G): ', s, ' alpha = ', alpha, 'gamma = ', gamma)
        print('edge density: ', fedge)
        return d, s, fedge
    else:
        d, s, fedge = rho_bipartite(G), edge_surplus_bipartite(G, alpha, gamma), fe_bipartite(G)
        print('n,m :', G.number_of_nodes(), G.number_of_edges())
        print('Bipartite graph statistics')
        print('rho(G): ', d)
        print('edge-surplus(G): ', s, ' alpha = ', alpha, 'gamma = ', gamma)
        print('edge density: ', fedge)
        return d, s, fedge

