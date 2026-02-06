from bioservices import KEGG

# Ignore warnings
def ignore_warnings():
    warnings.filterwarnings("ignore")

    
# Set random seed
def set_random_seed(seed=0):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


# Argument parser
def parse_arguments(labsm=0.3, wd=0.1, lr=0.0005,
                   n_hid=104, nheads=8, nlayers=3, 
                   cell_size=30, neighbor=20, egrn=True, 
                   output_file=None):
    """
    Parses command-line arguments for training a Graph Neural Network (GNN) on a gene-cell graph.
    
    Arguments:
        labsm (float): The rate of label smoothing to apply during training. Default is 0.3.
        wd (float): The weight decay parameter for regularization. Default is 0.1.
        lr (float): The learning rate for the optimizer. Default is 0.0005.
        n_hid (int): The number of nodes (or units) in the hidden layer. Default is 104.
        nheads (int): The number of attention heads used in multi-head attention layers. Default is 8.
        nlayers (int): The number of graph convolution layers in the GNN. Default is 3.
        cell_size (int): The feature dimension or size of a single cell in the gene-cell graph. Default is 30.
        neighbor (int): The number of neighboring nodes to consider for each node during graph convolution. Default is 20.
        egrn (bool): A flag to enable or disable Gene Regulatory Network (EGRN) functionality. Default is True.
        output_file (str): The path where the results will be saved. If not provided, the current directory will be used.
        
    Returns:
        tuple: A tuple containing the parsed arguments. These values will be used for training the model.
    """

    # If output_file is not provided, use the current working directory
    if output_file is None:
        output_file = os.getcwd() + '/Result/'

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description='Training GNN on gene cell graph')

    # Argument for label smoothing rate
    parser.add_argument('--labsm', type=float, default=labsm, 
                        help='The rate of Label Smoothing for reducing overfitting. Default is 0.3.')
    
    # Argument for weight decay (regularization)
    parser.add_argument('--wd', type=float, default=wd, 
                        help='The weight decay (L2 regularization) for the optimizer. Default is 0.1.')
    
    # Argument for learning rate
    parser.add_argument('--lr', type=float, default=lr, 
                        help='The learning rate for the optimizer. Default is 0.0005.')
    
    # Argument for the number of hidden units (nodes) in the hidden layers
    parser.add_argument('--n_hid', type=int, default=n_hid, 
                        help='The number of nodes (units) in the hidden layer. Default is 104.')
    
    # Argument for the number of attention heads in multi-head attention layers
    parser.add_argument('--nheads', type=int, default=nheads, 
                        help='The number of attention heads in multi-head attention layers. Default is 8.')
    
    # Argument for the number of graph convolution layers
    parser.add_argument('--nlayers', type=int, default=nlayers, 
                        help='The number of graph convolution layers in the model. Default is 3.')
    
    # Argument for the size (dimension) of a single cell in the graph
    parser.add_argument('--cell_size', type=int, default=cell_size, 
                        help='The feature dimension or size of a single cell. Default is 30.')
    
    # Argument for the number of neighboring nodes considered in graph convolution
    parser.add_argument('--neighbor', type=int, default=neighbor, 
                        help='The number of neighboring nodes considered for each node during graph convolution. Default is 20.')
    
    # Argument for enabling or disabling EGRN (Gene Regulatory Network)
    parser.add_argument('--egrn', type=bool, default=egrn, 
                        help='Enable or disable Gene Regulatory Network (EGRN) functionality. Default is True.')
    
    # Argument for output file path where results will be saved
    parser.add_argument('--output_file', type=str, default=output_file, 
                        help='Directory path where the results will be saved. Default is the current working directory + /Result/.')

    # Parse the arguments from the command line
    args = parser.parse_args()
    
    return (args.output_file, args.labsm, args.lr, args.wd, args.n_hid, 
            args.nheads, args.nlayers, args.cell_size, args.neighbor, args.egrn)


# Check if the folder exists, and create it if it does not
def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.makedirs(file_path)


def get_cancer_metastasis_genes():
    """
    Retrieve a list of gene symbols associated with cancer metastasis pathways 
    from the KEGG database.

    This function queries the KEGG database for a predefined set of cancer metastasis 
    related pathways, retrieves the associated genes for each pathway, and returns 
    them in a dictionary format.

    The cancer metastasis pathways considered are:
        - VEGF signaling pathway (hsa04370)
        - Focal adhesion (hsa04510)
        - Pathways in cancer (hsa05200)
        - PI3K-Akt signaling pathway (hsa04151)
        - MAPK signaling pathway (hsa04010)
        - Hippo signaling pathway (hsa04390)
        - ECM-receptor interaction (hsa04512)

    Returns:
        dict: A dictionary where keys are the pathway names (str) and values are 
              lists of gene symbols (list of str) associated with the respective 
              pathway.
              
              Example:
              {
                  'VEGF signaling pathway': ['Gene A', 'Gene B', 'Gene C', ...],
                  'Focal adhesion': ['Gene A', 'Gene B', 'Gene C', ...],
                  ...
              }

    """
    kegg = KEGG()

    # Define cancer metastasis pathways
    cancer_metastasis_pathways = {
        'VEGF signaling pathway': 'hsa04370',
        'Focal adhesion': 'hsa04510',
        'Pathways in cancer': 'hsa05200',
        'PI3K-Akt signaling pathway': 'hsa04151',
        'MAPK signaling pathway': 'hsa04010',
        "Hippo signaling pathway": "hsa04390",
        "ECM-receptor interaction": "hsa04512"
    }

    pathway_genes = {}

    # Iterate over each pathway
    for pathway_name, pathway_id in cancer_metastasis_pathways.items():
        # Retrieve pathway information
        pathway_info = kegg.get(pathway_id)
        parsed_pathway = kegg.parse(pathway_info)

        # Get and store the gene list
        genes = parsed_pathway['GENE']
        gene_symbols = []
        for gene_id, gene_info in genes.items():
            # Get the gene symbol
            gene_symbol = gene_info.split(' ')[0].split(';')[0]
            gene_symbols.append(gene_symbol)

        pathway_genes[pathway_name] = gene_symbols

    return pathway_genes
    

def initial_clustering(RNA_matrix, custom_n_neighbors=None, n_pcs=40, custom_resolution=None, use_rep=None, random_seed=0):
    """
    Perform clustering (Leiden clustering) on a gene-cell graph based on the input RNA expression matrix.

    This function takes a gene expression matrix, normalizes it, and performs Leiden clustering based on 
    the neighborhood graph of cells. The resolution and number of neighbors can be customized or automatically 
    determined based on the number of cells in the dataset.

    Parameters:
        RNA_matrix (numpy.ndarray or pandas.DataFrame): 
            The input gene expression matrix where rows represent genes and columns represent cells. 
            This is a **required parameter**.

        custom_n_neighbors (int, optional):
            The number of neighbors to consider when constructing the neighborhood graph. 
            If not provided, it will be automatically calculated based on the number of cells.
            Default is `None`, which triggers automatic calculation.

        n_pcs (int, optional):
            The number of principal components to use for neighborhood graph construction if `use_rep` is not provided. 
            Default is 40.

        custom_resolution (float, optional):
            The resolution for the Leiden clustering algorithm. If not provided, it will be automatically calculated based on the number of cells.
            Default is `None`, which triggers automatic calculation.

        use_rep (str or numpy.ndarray, optional):
            If provided, this will be used as a precomputed representation (embedding) for the neighborhood calculation. 
            If not provided, `n_pcs` will be used to compute the neighborhood graph.
            Default is `None`.

        random_seed (int, optional):
            The seed for the random number generator to ensure reproducibility. Default is 0.

    Returns:
        pandas.Series:
            A pandas Series containing the Leiden clustering labels for each cell.
            These labels indicate the cluster assignments for the cells in the input matrix.

    Notes:
        - `RNA_matrix` is a required parameter and must be provided by the user.
        - If `custom_n_neighbors` and `custom_resolution` are not provided, they will be calculated automatically 
          based on the number of cells in the dataset using the `segment_function`.
        - If `use_rep` is provided, it will be used for calculating the neighborhood graph instead of `n_pcs`.
    """
    print(
        '\tWhen the number of cells is less than or equal to 500, it is recommended to set the resolution value to 0.2.')
    print('\tWhen the number of cells is within the range of 500 to 5000, the resolution value should be set to 0.5.')
    print('\tWhen the number of cells is greater than 5000, the resolution value should be set to 0.8.')

    def segment_function(x):
        if x <= 500:
            return 0.2, 5
        elif x <= 5000:
            return 0.5, 10
        else:
            return 0.8, 15

    adata = ad.AnnData(RNA_matrix.transpose(), dtype='int32')

    # If the user did not provide a custom resolution or n_neighbors value, use the values calculated by segment_function
    if custom_resolution is None or custom_n_neighbors is None:
        resolution, n_neighbors = segment_function(adata.shape[0])
    else:
        resolution = custom_resolution
        n_neighbors = custom_n_neighbors

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Use the user-provided embedding if available, otherwise use n_pcs
    if use_rep is not None:
        adata.obsm['use_rep']=use_rep
        sc.pp.neighbors(adata, use_rep='use_rep', n_neighbors=n_neighbors)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    sc.tl.leiden(adata, resolution, random_state=random_seed)
    return adata.obs['leiden']
    
    
def subgraph_1(graph, seed, n_neighbors, node_sele_prob):
    total_matrix_size = 1 + np.cumprod(n_neighbors).sum()  # Number of nodes in the subgraph
    picked_nodes = {seed}  # One node in the batch
    last_layer_nodes = {seed}

    # Number of nodes selected in each layer. Initially, only the seed node is selected.
    to_pick = 1
    for n_neighbors_current in n_neighbors:  # Current layer neighbors
        to_pick = to_pick * n_neighbors_current
        neighbors = graph[list(last_layer_nodes), :].nonzero()[1]  # Find neighbors of last_layer_nodes

        neighbors_prob = node_sele_prob[list(neighbors)]
        neighbors = list(set(neighbors))  # Make all nodes from the last layer part of the neighbors set
        n_neigbors_real = min(
            to_pick,
            len(neighbors))  # Handle the case where the required number of neighbors is less than the actual number of neighbors
        if len(neighbors_prob) == 0:
            continue
        last_layer_nodes = set(
            np.random.choice(neighbors, n_neigbors_real, replace=False,
                             p=softmax(neighbors_prob)))  # Select non-repeated nodes from neighbors
        picked_nodes |= last_layer_nodes  # Update picked_nodes as last_layer_nodes ∪ picked_nodes
    indices = list(sorted(picked_nodes - {seed}))
    return indices


def softmax(x):
    return (np.exp(x) / np.exp(x).sum())

def subgraph_2(graph, seed, n_neighbors):
    picked_nodes = {seed}
    last_layer_nodes = {seed}

    for n_neighbors_current in n_neighbors:
        neighbors = graph[list(last_layer_nodes), :].nonzero()[1]
        neighbors = list(set(neighbors))
        n_neigbors_real = min(n_neighbors_current, len(neighbors))

        if len(neighbors) == 0:
            continue

        last_layer_nodes = set(np.random.choice(neighbors, n_neigbors_real, replace=False))
        picked_nodes |= last_layer_nodes

    indices = list(sorted(picked_nodes - {seed}))
    return indices

def batch_process_1(args):
    i, node_ids, RNA_matrix, neighbor, cell_size = args
    gene_indices_all = []
    dic = {}

    start_index = i * cell_size
    end_index = min(start_index + cell_size, len(node_ids))

    for node_index in node_ids[start_index:end_index]:
        rna_ = RNA_matrix[:, node_index].todense()
        rna_[rna_ < 5] = 0
        gene_indices = subgraph_1(RNA_matrix.transpose(), node_index, neighbor, np.squeeze(np.array(np.log(rna_ + 1))))
        dic[node_index] = {'g': gene_indices}
        gene_indices_all.extend(gene_indices)

    gene_indices_all = list(set(gene_indices_all))
    node_indices_all = node_ids[start_index:end_index]
    h = {'gene_index': gene_indices_all, 'cell_index': node_indices_all}
    return (h, dic)

def batch_process_2(args):
    i, node_ids, RNA_matrix, neighbor, cell_size = args
    gene_indices_all = []
    dic = {}

    start_index = i * cell_size
    end_index = min(start_index + cell_size, len(node_ids))

    for node_index in node_ids[start_index:end_index]:
        gene_indices = subgraph_2(RNA_matrix.transpose(), node_index, neighbor)
        dic[node_index] = {'g': gene_indices}
        gene_indices_all.extend(gene_indices)

    gene_indices_all = list(set(gene_indices_all))
    node_indices_all = node_ids[start_index:end_index]
    h = {'gene_index': gene_indices_all, 'cell_index': node_indices_all}
    return (h, dic)


def batch_select_whole(RNA_matrix, label, neighbor=[20], cell_size=30):
    print('Partitioning data into batches based on sample type.')
    dic = {}

    # Randomly shuffle cell indices
    shuffled_indices = np.random.choice(len(label), size=len(label), replace=False)

    # Create a mapping from shuffled index to original index
    index_mapping = {new_idx: original_idx for original_idx, new_idx in enumerate(shuffled_indices)}

    # Get shuffled labels using the shuffled index mapping
    shuffled_labels = [label[index_mapping[i]] for i in range(len(label))]

    # Separate samples into type P and M
    p_node_ids = [index_mapping[i] for i, l in enumerate(shuffled_labels) if l == 'P']
    m_node_ids = [index_mapping[i] for i, l in enumerate(shuffled_labels) if l == 'M']

    # Calculate the number of batches needed for each sample type
    n_batch_p = math.ceil(len(p_node_ids) / cell_size)
    n_batch_m = math.ceil(len(m_node_ids) / cell_size)

    with mp.Pool(processes=48) as pool:
        # Process P samples using batch_process_2
        tasks_p = [(i, p_node_ids, RNA_matrix, neighbor, cell_size) for i in range(n_batch_p)]
        results_p = list(tqdm(pool.imap_unordered(batch_process_2, tasks_p), total=n_batch_p, desc="Processing P samples"))

        # Process M samples using batch_process_1
        tasks_m = [(i, m_node_ids, RNA_matrix, neighbor, cell_size) for i in range(n_batch_m)]
        results_m = list(tqdm(pool.imap_unordered(batch_process_1, tasks_m), total=n_batch_m, desc="Processing M samples"))

    # Combine results
    results = results_p + results_m
    indices_ss = [res[0] for res in results]
    for res in results:
        dic.update(res[1])

    all_cell_indices = [index for batch in indices_ss for index in batch['cell_index']]
    
    # Return the node indices in original order
    return indices_ss, all_cell_indices, dic
