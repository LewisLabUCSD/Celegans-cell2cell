import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def get_element_num(dataframe):
    lr_elements = 0
    for i in dataframe['Ligand']:
        lr_elements += 1
    return lr_elements


def get_coordinates(radius, G, sorting_feature=None):   
    elements = list(set(G.nodes()))
    n_elements = len(elements)
    
    if sorting_feature is not None:
        sorting_list = [G.nodes[n][sorting_feature] for n in elements]
        elements = [x for _, x in sorted(zip(sorting_list, elements))]
    
    coordinates = dict()
    
    for i, node in enumerate(elements):
        theta = (2*np.pi / n_elements) * i
        x = radius*np.cos(theta)
        y = radius*np.sin(theta)
        coordinates[node] = (x, y)
    return coordinates


def get_label_coordinates(radius, G, sorting_feature=None):   
    elements = list(set(G.nodes()))
    n_elements = len(elements)
    
    if sorting_feature is not None:
        sorting_list = [G.nodes[n][sorting_feature] for n in elements]
        elements = [x for _, x in sorted(zip(sorting_list, elements))]
    
    coordinates = dict()
    
    for i, node in enumerate(elements):
        theta = (2*np.pi / n_elements) * i
        theta2 = (360 / n_elements) * i
        x = radius*np.cos(theta)
        y = radius*np.sin(theta)
        coordinates[node] = (x, y, theta2)
    return coordinates


def determine_small_radius(coordinate_dict):
    keys = list(coordinate_dict.keys())
    circle1 = np.asarray(coordinate_dict[keys[0]])
    circle2 = np.asarray(coordinate_dict[keys[1]])
    diff = circle1 - circle2
    distance = np.sqrt(diff[0]**2 + diff[1]**2)
    return distance / 2.0


def get_node_colors(G, coloring_feature=None, cmap='viridis'):
    if coloring_feature is None:
        raise ValueError('Node feature not specified!')
    
    # Get what features we have in the network G
    features = set()
    for n in G.nodes():
        features.add(G.nodes[n][coloring_feature])
    
    # Generate colors for each feature
    NUM_COLORS = len(features)
    cm = plt.get_cmap(cmap)
    feature_colors = dict()
    for i, f in enumerate(features):
        feature_colors[f] = cm(1.*i/NUM_COLORS)  # color will now be an RGBA tuple    
    
    # Map feature colors into nodes
    node_colors = dict()
    for n in G.nodes():
        feature = G.nodes[n][coloring_feature]
        node_colors[n] = feature_colors[feature]
    return node_colors, feature_colors
    
    
def get_color_edges_continuous(G, coloring_feature='weight', min_val=0, max_val=None, cmap='tab20'):
    if coloring_feature is None:
        raise ValueError('Node feature not specified!')
    #Get the num of unique features in coloring_feature:
    features = []
    for f in G.edges():
        features.append(G.edges[f][coloring_feature])
        
    if max_val is None:
        max_val = max(features)
    norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val, clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    
    edge_colors = dict()
    for e in G.edges():
        feature = G.edges[e][coloring_feature]

        edge_colors[e] = mapper.to_rgba(feature)

    steps = max_val/10.
    feature_colors = {v : mapper.to_rgba(v) for v in np.arange(0, max_val + steps, steps)}
    return edge_colors, feature_colors



        
    