import networkx as nx
import CreateSpecies as cs


def connectivity_index_analysis():
    ferret_nasal_ds = cs.create_ferret_nasal(is_connectivity_index=True)
    A = nx.nx_agraph.to_agraph(ferret_nasal_ds.max_connectivity_index_graph)
    A.layout()
    A.draw('/Users/maxbagga/Box/Mallard Phylogeny Data/Ferret Nasal Results/Ferret Nasal Max Connectivity Index Graph.png')
    print(ferret_nasal_ds.average_connectivity_index)

    guinea_pig_nasal_ds = cs.create_guinea_pig_nasal(is_connectivity_index=True)
    A = nx.nx_agraph.to_agraph(guinea_pig_nasal_ds.max_connectivity_index_graph)
    A.layout()
    A.draw('/Users/maxbagga/Box/Mallard Phylogeny Data/Guinea Pig Nasal Results/Guinea Pig Nasal Max Connectivity Index '
           'Graph.png')
    print(guinea_pig_nasal_ds.average_connectivity_index)

    pig_nasal_ds = cs.create_pig_nasal(is_connectivity_index=True)
    A = nx.nx_agraph.to_agraph(pig_nasal_ds.max_connectivity_index_graph)
    A.layout()
    A.draw('/Users/maxbagga/Box/Mallard Phylogeny Data/Pig Nasal Results/Pig Nasal Max Connectivity Index graph.png')
    print(pig_nasal_ds.average_connectivity_index)

    quail_tracheal_ds = cs.create_quail_tracheal(is_connectivity_index=True)
    A = nx.nx_agraph.to_agraph(quail_tracheal_ds.max_connectivity_index_graph)
    A.layout()
    A.draw('/Users/maxbagga/Box/Mallard Phylogeny Data/Quail Tracheal Results/Quail Tracheal Max Connectivity Index '
           'Graph.png')
    print(quail_tracheal_ds.average_connectivity_index)