from .fibheap import FibonacciHeap as FibHeap
import networkx as nx
import copy

# class DSFuncs(object):
    # def __init__(self, arg):
    #   super(DenseSubgraphFuncs, self).__init__()
    #   self.arg = arg
    # @staticmethod
def check_undirected_graph(G):
    if isinstance(G, nx.Graph):
        return

    raise NotImplementedError('Only accept networkx undirected graph (nx.Graph).')

# @staticmethod
def exact_densest(hyper_list: list):
    """
    Goldberg's exact max flow algorithm.

    Parameters
    ----------
    G: undirected, graph (networkx).

    Returns
    -------
    Sstar: list, subset of nodes corresponding to densest subgraph.
    opt: float, density of Sstar induced subgraph.

    """
    # self.check_undirected_graph(G)

    if isinstance(hyper_list, nx.Graph):
        return exact_densest_from_graph(hyper_list)

    weight = len(hyper_list[0])
    node_w = dict()
    for e in hyper_list:
        for node in e:
            node_w[node] = 1.0 if node not in node_w else node_w[node]+1

    inf = len(hyper_list)*10.0

    m = float(len(hyper_list))
    n = float(len(node_w))
    minD = m / n  # rho^* >= m/n since V is a feasible solution
    maxD = m
    # a tighter upper bound
    # degree_seq = [d for n,d in G.degree()]
    # maxD = (max(degree_seq)-1 )/2

    opt = minD
    Sstar = node_w.keys()
    # print(maxD, minD)
    if minD == maxD:
        return Sstar, maxD
    while maxD - minD > 1 / n ** 2:  # binary search
        query = (maxD + minD) / 2
        # print('Query value is ',query )
        H = create_flow_network_from_list(hyper_list, query, node_w, inf, weight)
        solution = nx.minimum_cut(H, 's', 't', capacity='capacity')  # unspecified behavior
        
        cut = solution[1][0]
        # print(solution[0], cut)
        if cut == {'s'}:
            maxD = query  # this means there is no subgraph S such that the degree density is at least query
        else:
            #             print('Found denser subgraph!')
            minD = query
            Sstar = cut
            opt = query
    Sstar = list(set(Sstar)&set(node_w.keys()))
    return Sstar, opt

# @staticmethod
def exact_densest_from_graph(G: nx.Graph):
    """
    Goldberg's exact max flow algorithm.

    Parameters
    ----------
    G: undirected, graph (networkx).

    Returns
    -------
    Sstar: list, subset of nodes corresponding to densest subgraph.
    opt: float, density of Sstar induced subgraph.

    """
    check_undirected_graph(G)

    m = G.number_of_edges()
    n = G.number_of_nodes()
    minD = m / n  # rho^* >= m/n since V is a feasible solution
    maxD = (n - 1) / 2
    # a tighter upper bound
    # degree_seq = [d for n,d in G.degree()]
    # maxD = (max(degree_seq)-1 )/2

    opt = 0
    Sstar = G.nodes()
    if minD == maxD:
        return Sstar, maxD
    while maxD - minD > 1 / n ** 2:  # binary search
        query = (maxD + minD) / 2
        # print('Query value is ',query )
        H = create_flow_network_from_graph(G, query)
        solution = nx.minimum_cut(H, 's', 't', capacity='capacity')  # unspecified behavior
        #         print(solution[0])
        cut = solution[1][0]
        #         print(cut)
        if cut == {'s'}:
            maxD = query  # this means there is no subgraph S such that the degree density is at least query
        else:
            #             print('Found denser subgraph!')
            minD = query
            Sstar = cut
            opt = query
    Sstar = list(set(Sstar)&set(G.nodes()))
    return Sstar, opt


# @staticmethod
def greedy_charikar(Ginput, weight=None):
    """
    Charikar's 1/2 greedy algorithm

    This function greedily removes nodes, and ouputs a set of nodes that optimizes
    the objective rho(S) = e[S]/|S|


    Parameters
    ----------
    Ginput: undirected, graph (networkx)
    weight: str that specify the edge attribute name of edge weight; None if the graph is unweighted

    Returns
    -------
    H: list, subset of nodes  corresponding to densest subgraph
    rho: float, density of H induced subgraph

    """

    H, rho, _ = flowless(Ginput, 1, weight=weight)
    return H, rho


# @staticmethod
def flowless(hyper_list, T, weight=None):
    """
    Implementation of greedy++ algorithm proposed in [1]. It efficiently compute almost densest subgraph without max
    flow.

    Parameters
    ----------
    G: undirected, graph (networkx).
    T: int, number of iterations.
    weight: str that specify the edge attribute name of edge weight; None if the graph is unweighted.

    Returns
    ----------
    Hstar: list, subset of nodes corresponding to densest subgraph.
    rhostar: float, density of Hstar induced subgraph.
    loadsstar: dict, loads for nodes.

     References
     ----------
     [1] Digvijay Boob, Yu Gao, Richard Peng, Saurabh Sawlani, Charalampos Tsourakakis,
     Di Wang, and Junxing Wang. 2020. Flowless: Extracting Densest Subgraphs Without Flow Computations. In
     Proceedings of The Web Conference 2020 (WWW '20).
    """
    # self.check_undirected_graph(G)
    if isinstance(hyper_list, nx.Graph):
        return flowless_from_graph(hyper_list, T, weight)

    loads = dict()
    loadsstar = loads
    # H = copy.deepcopy(G)
    node_dict, fibheap, total_degree = init_heap_flowless(hyper_list, loads, weight=weight)
    rhostar = float(len(hyper_list)) / len(node_dict)
    Hstar = list(node_dict.keys())
    for i in range(T):
        if i>0:
            node_dict, fibheap, total_degree = init_heap_flowless(hyper_list, l, weight=weight)
        H, tmp, l = greedy_helper(hyper_list, node_dict, fibheap, total_degree, weight=weight)
        print('iteration ' + repr(i + 1))
        if tmp > rhostar:
            rhostar = tmp
            Hstar = H
            loadsstar = l
    return Hstar, rhostar, loadsstar


# @staticmethod
def flowless_from_graph(G, T, weight=None):
    """
    Implementation of greedy++ algorithm proposed in [1]. It efficiently compute almost densest subgraph without max
    flow.

    Parameters
    ----------
    G: undirected, graph (networkx).
    T: int, number of iterations.
    weight: str that specify the edge attribute name of edge weight; None if the graph is unweighted.

    Returns
    ----------
    Hstar: list, subset of nodes corresponding to densest subgraph.
    rhostar: float, density of Hstar induced subgraph.
    loadsstar: dict, loads for nodes.

     References
     ----------
     [1] Digvijay Boob, Yu Gao, Richard Peng, Saurabh Sawlani, Charalampos Tsourakakis,
     Di Wang, and Junxing Wang. 2020. Flowless: Extracting Densest Subgraphs Without Flow Computations. In
     Proceedings of The Web Conference 2020 (WWW '20).
    """
    # self.check_undirected_graph(G)

    loads = dict()
    loadsstar = loads
    H = copy.deepcopy(G)
    rhostar = H.number_of_edges() / H.number_of_nodes()
    Hstar = H
    for i in range(T):
        node_dict, fibheap, total_degree = init_heap_flowless_from_graph(G, loads, weight=weight)
        H, tmp, loads = greedy_helper_from_graph(G, node_dict, fibheap, total_degree, weight=weight)
        # print('iteration ' + repr(i + 1))
        if tmp > rhostar:
            rhostar = tmp
            Hstar = H
            loadsstar = loads
    return Hstar, rhostar, loadsstar




def create_flow_network_from_graph(G, query):
    m = G.number_of_edges()
    G = nx.DiGraph(G)
    H = G.copy()
    H.add_node('s')
    H.add_node('t')
    for e in G.edges():
        H.add_edge(e[0], e[1], capacity=1)

    for v in G.nodes():
        H.add_edge('s', v, capacity=m)
        H.add_edge(v, 't', capacity=m + 2 * query - G.in_degree(v))
    return H


def create_flow_network_from_list(hyper_list, query, node_w, inf, weight):
    


    # m = G.number_of_edges()
    # G = nx.DiGraph(G)
    H = nx.DiGraph()
    H.add_node('s')
    H.add_node('t')
    for e in hyper_list:
        for node in e:
            H.add_edge(node, tuple(e), capacity=1.0/weight)
            H.add_edge(tuple(e), node, capacity=inf)

    for v in node_w:
        H.add_edge('s', v, capacity=node_w[v]/weight)
        H.add_edge(v, 't', capacity=query)
    return H

def init_heap_flowless(hyper_list, loads=None, weight=None):
    """
    Parameters
    ----------
    G: undirected, graph (networkx)
    loads: dict, initial load for each node. Used in algorithm flowless, see
    weight: str that specify the edge attribute name of edge weight; None if the graph is unweighted

    Returns
    ----------
    node_dict: dict, node id as key, tuple (neighbor list, heap node) as value. Here heap node is a
    pointer to the corresponding node in fibheap.
    fibheap: FibonacciHeap, support fast extraction of min degree node and value change.
    total_weight: edge weight sum.
    """
    node_dict = dict()
    fibheap = FibHeap()
    total_degree = 0

    for e_index in range(len(hyper_list)):
        e = hyper_list[e_index]
        for node in e:
            if node not in node_dict:
                node_dict[node] = (list(), fibheap.insert(0, node))
            fibheap.decrease_key(node_dict[node][1], node_dict[node][1].key + 1)
            node_dict[node][0].append(e_index)
    # for node in G.nodes:
    #     node_dict[node] = (list(), fibheap.insert(0, node))
    #     for neighbor in G[node]:
    #         if weight is None:
    #             edge_w = 1
    #         else:
    #             edge_w = G[node][neighbor][weight]
    #         fibheap.decrease_key(node_dict[node][1], node_dict[node][1].key + edge_w)
    #         node_dict[node][0].append(neighbor)
    #         total_degree += edge_w
    total_weight = len(hyper_list)

    for node in loads:
        fibheap.decrease_key(node_dict[node][1], node_dict[node][1].key + loads[node])

    return node_dict, fibheap, total_weight


def init_heap_flowless_from_graph(G, loads=None, weight=None):
    """
    Parameters
    ----------
    G: undirected, graph (networkx)
    loads: dict, initial load for each node. Used in algorithm flowless, see
    weight: str that specify the edge attribute name of edge weight; None if the graph is unweighted

    Returns
    ----------
    node_dict: dict, node id as key, tuple (neighbor list, heap node) as value. Here heap node is a
    pointer to the corresponding node in fibheap.
    fibheap: FibonacciHeap, support fast extraction of min degree node and value change.
    total_weight: edge weight sum.
    """
    node_dict = dict()
    fibheap = FibHeap()
    total_degree = 0
    for node in G.nodes():
        node_dict[node] = (list(), fibheap.insert(0, node))
        for neighbor in G[node]:
            if weight is None:
                edge_w = 1
            else:
                edge_w = G[node][neighbor][weight]
            fibheap.decrease_key(node_dict[node][1], node_dict[node][1].key + edge_w)
            node_dict[node][0].append(neighbor)
            total_degree += edge_w
    total_weight = total_degree / 2

    for node in loads:
        fibheap.decrease_key(node_dict[node][1], node_dict[node][1].key + loads[node])

    return node_dict, fibheap, total_weight


def greedy_helper_from_graph(G, node_dict, fib_heap, total_degree, weight=None):
    """
    Greedy peeling algorithm. Peel nodes iteratively based on their current degree.

    Parameters
    ----------
    G: undirected, graph (networkx)
    node_dict: dict, node id as key, tuple (neighbor list, heap node) as value. Here heap node is a
    pointer to the corresponding node in fibheap.
    fibheap: FibonacciHeap, support fast extraction of min degree node and value change.
    total_weight: edge weight sum.
    weight: str that specify the edge attribute name of edge weight; None if the graph is unweighted.

    Returns
    ----------
    H: list, subset of nodes corresponding to densest subgraph.
    max_avg: float, density of H induced subgraph.
    new_loads: dict, new loads for nodes, only used for the flowless algorithm when T>1.
    """
    n = G.number_of_nodes()
    avg_degree = total_degree / n
    H = list(G)
    max_avg = avg_degree
    new_loads = dict()

    for i in range(n - 1):
        # find min node from graph (remove from heap)
        to_remove = fib_heap.extract_min()
        node_to_remove = to_remove.value
        degree_to_remove = to_remove.key
        new_loads[node_to_remove] = degree_to_remove

        # for every neighbor node this min node have
        for neighbor in node_dict[node_to_remove][0]:
            edge_w = 1 if weight is None else G[node_to_remove][neighbor][weight]

            # here the key can be actually increased
            if neighbor != node_to_remove:
                fib_heap.decrease_key(node_dict[neighbor][1], node_dict[neighbor][1].key - edge_w)
                node_dict[neighbor][0].remove(node_to_remove)
            total_degree -= edge_w

        del node_dict[node_to_remove]
        avg_degree = total_degree / (n - i - 1)
        if max_avg < avg_degree:
            max_avg = avg_degree
            H = list(node_dict.keys())

    return H, max_avg, new_loads

def greedy_helper(hyper_list, node_dict, fib_heap, total_weight, weight=None):
    """
    Greedy peeling algorithm. Peel nodes iteratively based on their current degree.

    Parameters
    ----------
    G: undirected, graph (networkx)
    node_dict: dict, node id as key, tuple (neighbor list, heap node) as value. Here heap node is a
    pointer to the corresponding node in fibheap.
    fibheap: FibonacciHeap, support fast extraction of min degree node and value change.
    total_weight: edge weight sum.
    weight: str that specify the edge attribute name of edge weight; None if the graph is unweighted.

    Returns
    ----------
    H: list, subset of nodes corresponding to densest subgraph.
    max_avg: float, density of H induced subgraph.
    new_loads: dict, new loads for nodes, only used for the flowless algorithm when T>1.
    """
    n = len(node_dict.keys())
    avg_degree = total_weight / n
    H = list(node_dict.keys())
    max_avg = avg_degree
    new_loads = dict()

    for i in range(n - 1):
        # find min node from graph (remove from heap)
        to_remove = fib_heap.extract_min()
        node_to_remove = to_remove.value
        degree_to_remove = to_remove.key
        new_loads[node_to_remove] = degree_to_remove

        for e_index in node_dict[node_to_remove][0]:
            e = hyper_list[e_index]
            for neighbor in e:
                if neighbor != node_to_remove:
                    fib_heap.decrease_key(node_dict[neighbor][1], node_dict[neighbor][1].key - 1)
                    node_dict[neighbor][0].remove(e_index)
            total_weight -= 1

        del node_dict[node_to_remove]
        avg_degree = total_weight / (n - i - 1)
        if max_avg < avg_degree:
            max_avg = avg_degree
            H = list(node_dict.keys())

    return H, max_avg, new_loads
