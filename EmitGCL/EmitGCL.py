import os
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
from tqdm import tqdm
from torch.utils.data import SubsetRandomSampler
from loss_function import *
from conv import *
from scipy.sparse import csc_matrix, csr_matrix
from scipy.io import mmwrite,mmread


class GNN_from_raw(nn.Module):
    def __init__(self, in_dim, n_hid, num_types, num_relations, n_heads, n_layers, dropout=0.2, conv_name='hgt',
                 prev_norm=True, last_norm=True):
        super(GNN_from_raw, self).__init__()
        self.gcs = nn.ModuleList()
        self.num_types = num_types
        self.in_dim = in_dim
        self.n_hid = n_hid
        self.adapt_ws = nn.ModuleList()
        self.drop = nn.Dropout(dropout)
        self.embedding1 = nn.ModuleList()

        # Initialize MLP weight matrices
        for ti in range(num_types):
            self.embedding1.append(nn.Linear(in_dim[ti], 256))

        for t in range(num_types):
            self.adapt_ws.append(nn.Linear(256, n_hid))

        # Initialize graph convolution layers
        for l in range(n_layers - 1):
            self.gcs.append(
                GeneralConv(conv_name, n_hid, n_hid, num_types, num_relations, n_heads, dropout, use_norm=prev_norm))
        self.gcs.append(
            GeneralConv(conv_name, n_hid, n_hid, num_types, num_relations, n_heads, dropout, use_norm=last_norm))

    def encode(self, x, t_id):
        h1 = F.relu(self.embedding1[t_id](x))
        return h1

    def forward(self, node_feature, node_type, edge_index, edge_type):
        node_embedding = []
        for t_id in range(self.num_types):
            node_embedding += list(self.encode(node_feature[t_id], t_id))

        node_embedding = torch.stack(node_embedding)
        # Initialize result matrix
        res = torch.zeros(node_embedding.size(0), self.n_hid).to(node_feature[0].device)

        # Process each node type
        for t_id in range(self.num_types):
            idx = (node_type == int(t_id))
            if idx.sum() == 0:
                continue
            # Update result matrix
            res[idx] = torch.tanh(self.adapt_ws[t_id](node_embedding[idx]))

        # Apply dropout to the result matrix
        meta_xs = self.drop(res)
        del res

        # Initialize a container for the final layer's attention weights
        self.final_layer_attention_weights = {}
        # Iterate through graph convolution layers
        for i, gc in enumerate(self.gcs):
            meta_xs = gc(meta_xs, node_type, edge_index, edge_type)
            # Only store the attention weights from the final layer
            if i == len(self.gcs) - 1:
                self.final_layer_attention_weights = gc.edge_attention_weights

        return meta_xs


class Net(nn.Module):
    def __init__(self, dim_in, dim_out):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(dim_in, dim_out)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        return x
        
        
class NodeDimensionReduction(nn.Module):
    def __init__(self, RNA_matrix, indices, ini_p1, n_hid, n_heads, n_layers, labsm, lr, 
                 wd, device, num_types=2, num_relations=1, epochs=1):
        super(NodeDimensionReduction, self).__init__()
        self.RNA_matrix = RNA_matrix
        self.indices = indices
        self.ini_p1 = ini_p1
        self.in_dim = [RNA_matrix.shape[0], RNA_matrix.shape[1]]
        self.n_hid = n_hid
        self.num_types = num_types
        self.num_relations = num_relations
        self.n_heads = n_heads
        self.n_layers = n_layers
        self.labsm = labsm
        self.lr = lr
        self.wd = wd
        self.device = device
        self.epochs = epochs
        self.LabSm = LabelSmoothing(self.labsm)

        self.gnn = GNN_from_raw(in_dim=self.in_dim,
                                n_hid=self.n_hid,
                                num_types=self.num_types,
                                num_relations=self.num_relations,
                                n_heads=self.n_heads,
                                n_layers=self.n_layers,
                                dropout=0.3).to(self.device)

        self.optimizer = torch.optim.AdamW(self.gnn.parameters(), lr=self.lr, weight_decay=self.wd)
        self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(self.optimizer, 'min', factor=0.5, patience=5,
                                                                    verbose=True)

    def train_model(self, n_batch):
        print('The training process for the NodeDimensionReduction model has started. Please wait.')
        self.subgraph_attention_weights = {}

        for epoch in tqdm(range(self.epochs)):
            total_loss_kl1 = 0
            total_loss_cluster = 0
            all_cell_emb = []

            for batch_id in np.arange(n_batch):
                gene_index = self.indices[batch_id]['gene_index']
                cell_index = self.indices[batch_id]['cell_index']
                gene_feature = self.RNA_matrix[list(gene_index),]
                cell_feature = self.RNA_matrix[:, list(cell_index)].T
                gene_feature = torch.tensor(np.array(gene_feature.todense()), dtype=torch.float32).to(self.device)
                cell_feature = torch.tensor(np.array(cell_feature.todense()), dtype=torch.float32).to(self.device)

                node_feature = [cell_feature, gene_feature]
                gene_cell_sub = self.RNA_matrix[list(gene_index),][:, list(cell_index)]
                gene_cell_edge_index1 = list(np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1]) + list(
                    np.nonzero(gene_cell_sub)[1])
                gene_cell_edge_index2 = list(np.nonzero(gene_cell_sub)[1]) + list(
                    np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1])
                gene_cell_edge_index = torch.LongTensor([gene_cell_edge_index1, gene_cell_edge_index2]).to(self.device)

                edge_index = gene_cell_edge_index
                node_type = torch.LongTensor(np.array(
                    list(np.zeros(len(cell_index))) + list(np.ones(len(gene_index))))).to(self.device)
                edge_type = torch.LongTensor(np.array(list(np.zeros(np.nonzero(gene_cell_sub)[0].shape[0])) + list(
                    np.ones(np.nonzero(gene_cell_sub)[1].shape[0])))).to(self.device)
                l = torch.LongTensor(np.array(self.ini_p1)[[cell_index]]).to(self.device)

                node_rep = self.gnn.forward(node_feature, node_type,
                                            edge_index,
                                            edge_type).to(self.device)
                
                self.subgraph_attention_weights[batch_id] = self.gnn.final_layer_attention_weights
                
                cell_emb = node_rep[node_type == 0]
                gene_emb = node_rep[node_type == 1]
                
                all_cell_emb.append(cell_emb)
                
                decoder1 = torch.mm(gene_emb, cell_emb.t())
                gene_cell_sub = torch.tensor(np.array(gene_cell_sub.todense()), dtype=torch.float32).to(self.device)

                logp_x1 = F.log_softmax(decoder1, dim=-1)
                p_y1 = F.softmax(gene_cell_sub, dim=-1)

                loss_kl1 = F.kl_div(logp_x1, p_y1, reduction='mean')
                loss_cluster = self.LabSm(cell_emb, l)

                batch_loss = loss_kl1 + loss_cluster

                self.optimizer.zero_grad()
                batch_loss.backward()
                self.optimizer.step()

                total_loss_kl1 += loss_kl1.item()
                total_loss_cluster += loss_cluster.item()

            total_loss_kl1 /= n_batch
            total_loss_cluster /= n_batch

            all_cell_emb = torch.cat(all_cell_emb, dim=0)
            all_cell_pre = all_cell_emb.argmax(dim=1)
            # Count the occurrences of each unique value in all_cell_pre
            unique_values, counts = torch.unique(all_cell_pre, return_counts=True)
            # Print the unique values and their counts
            # print("Unique values in all_cell_pre and their counts:")
            # for value, count in zip(unique_values, counts):
                # print(f"Value: {value.item()}, Count: {count.item()}")

            print(f'Epoch {epoch + 1}:')
            print(f'  KL Loss: {total_loss_kl1}')
            print(f'  Cluster Loss: {total_loss_cluster}')

        print('The training for the NodeDimensionReduction model has been completed.')
        return self.gnn
        
        
class EmitGCL(nn.Module):
    def __init__(self, gnn, labsm, n_hid, n_batch, device, lr, wd, pathway_genes, pathway_directions, gene_names, sample_type, sample_list, num_epochs=1):
        super(EmitGCL, self).__init__()
        self.lr = lr
        self.wd = wd
        self.gnn = gnn
        self.n_hid = n_hid
        self.n_batch = n_batch
        self.device = device
        self.num_epochs = num_epochs
        self.net = Net(2 * self.n_hid, self.n_hid).to(self.device)
        self.gnn_optimizer = torch.optim.AdamW(self.gnn.parameters(), lr=self.lr, weight_decay=self.wd)
        self.net_optimizer = torch.optim.AdamW(self.net.parameters(), lr=1e-2)
        self.labsm = labsm
        self.LabSm = LabelSmoothing(self.labsm)
        
        # Initialize PathwayUCellLoss
        gene_name_to_index = {gene: i for i, gene in enumerate(gene_names)}
        self.pathway_loss_fn = PathwayUCellLoss(pathway_genes, gene_name_to_index, pathway_directions, device)  
        
        # Store the sample list for sequential sampling
        self.sample_list = sample_list
        
    def forward(self, indices, RNA_matrix, ini_p1, sample_type, nodes_id):
        cluster_l = list()
        cluster_kl_l = list()
        sim_l = list()
        nmi_l = list()
        ini_p1 = np.array(ini_p1)
        
        self.subgraph_attention_weights = {}
        
        for epoch in range(self.num_epochs):
            total_loss_kl1 = 0
            total_loss_cluster = 0
            total_ucell_loss = 0
            total_contrastive_loss = torch.tensor(0.0, device=self.device)  # Ensure the initial value is a tensor
            valid_contrastive_loss_count = 0
            all_cell_emb = []
            
            for batch_id in tqdm(np.arange(self.n_batch)):
                gene_index = indices[batch_id]['gene_index']
                cell_index = indices[batch_id]['cell_index']
                gene_feature = RNA_matrix[list(gene_index),]
                cell_feature = RNA_matrix[:, list(cell_index)].T
                gene_feature = torch.tensor(np.array(gene_feature.todense()), dtype=torch.float32).to(self.device)
                cell_feature = torch.tensor(np.array(cell_feature.todense()), dtype=torch.float32).to(self.device)
                node_feature = [cell_feature, gene_feature]
                gene_cell_sub = RNA_matrix[list(gene_index),][:, list(cell_index)]
                gene_cell_edge_index1 = list(np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1]) + list(
                    np.nonzero(gene_cell_sub)[1])
                gene_cell_edge_index2 = list(np.nonzero(gene_cell_sub)[1]) + list(
                    np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1])
                gene_cell_edge_index = torch.LongTensor([gene_cell_edge_index1, gene_cell_edge_index2]).to(self.device)
                node_type = torch.LongTensor(np.array(
                    list(np.zeros(len(cell_index))) + list(np.ones(len(gene_index))))).to(self.device)
                edge_type = torch.LongTensor(np.array(list(np.zeros(np.nonzero(gene_cell_sub)[0].shape[0])) + list(
                    np.ones(np.nonzero(gene_cell_sub)[1].shape[0])))).to(self.device)
                l = torch.LongTensor(np.array(ini_p1)[[cell_index]]).to(self.device)
                node_rep = self.gnn.forward(node_feature, node_type,
                                            gene_cell_edge_index,
                                            edge_type).to(self.device)
                
                self.subgraph_attention_weights[batch_id] = self.gnn.final_layer_attention_weights
                
                cell_emb = node_rep[node_type == 0]
                gene_emb = node_rep[node_type == 1]

                all_cell_emb.append(cell_emb)

                decoder1 = torch.mm(gene_emb, cell_emb.t())
                gene_cell_sub = torch.tensor(np.array(gene_cell_sub.todense()), dtype=torch.float32).to(self.device)
                logp_x1 = F.log_softmax(decoder1, dim=-1)
                p_y1 = F.softmax(gene_cell_sub, dim=-1)

                loss_kl1 = F.kl_div(logp_x1, p_y1, reduction='mean')
                total_loss_kl1 += loss_kl1
                loss_cluster = self.LabSm(cell_emb, l)
                total_loss_cluster += loss_cluster

            all_cell_emb = torch.cat(all_cell_emb, dim=0)
            all_cell_pre = all_cell_emb.argmax(dim=1)
            
            # Sequentially calculate ucell_loss and contrastive loss for sample pairs
            for i in range(len(self.sample_list) - 1):
                sample_1 = self.sample_list[i]
                sample_2 = self.sample_list[i + 1]
                
                # Create masks for the two samples in the current pair
                mask_1 = (np.array(sample_type) == sample_1)
                mask_2 = (np.array(sample_type) == sample_2)
                
                # Combine the two masks
                sample_pair_mask = mask_1 | mask_2
                
                print("RNA_matrix shape:", RNA_matrix.shape)
                
                RNA_matrix_pair = RNA_matrix[:,nodes_id][:, sample_pair_mask]
                labels_pair = all_cell_pre[sample_pair_mask]
                sample_type_pair = np.array(sample_type)[sample_pair_mask]
                sample_pair_emb = all_cell_emb[sample_pair_mask]
                # Compute ucell_loss and indices for high and low scores
                ucell_loss_tensor, high_score_indices, low_score_indices = self.pathway_loss_fn(
                    RNA_matrix_pair.transpose(), labels_pair, sample_type_pair,sample_2
                )
                
                total_ucell_loss += ucell_loss_tensor
                
                # Compute contrastive loss for each label in the current sample pair
                for label in high_score_indices:
                    label_mask = (labels_pair == label)
#                     print("label_mask.shape:", label_mask.shape)
#                     print("sample_type_pair.shape:", sample_type_pair.shape)
                    reference_mask = (sample_type_pair == sample_1) & label_mask.cpu().numpy()
#                     print("sample_pair_emb.shape:", sample_pair_emb.shape)
#                     print("reference_mask.shape:", reference_mask.shape)

                    z_reference = sample_pair_emb[reference_mask]
                    
                    high_score_features = sample_pair_emb[label_mask][high_score_indices[label]]
                    low_score_features = sample_pair_emb[label_mask][low_score_indices[label]]
                    
#                     other_cells_mask = ~np.isin(np.arange(len(labels_pair)), np.concatenate([high_score_indices[label], low_score_indices[label]]))
#                     z_other = sample_pair_emb[other_cells_mask]

                    # Define the mask to select cells that are not of the current label and belong to sample_2
                    other_cells_mask = (labels_pair != label) & (sample_type_pair == sample_2)
                    # Select the embeddings for these "other" cells
                    z_other = sample_pair_emb[other_cells_mask]

                    reference_indices = np.where(reference_mask)[0]
                    other_cells_indices = np.where(other_cells_mask)[0]
                    overlapping_indices = np.intersect1d(reference_indices, other_cells_indices)
                    print(f"Overlapping indices: {overlapping_indices}")
                    
                    # Sample other cells
                    sample_size = int(len(z_other) * 1.0)
                    sampler = torch.utils.data.SubsetRandomSampler(torch.randperm(len(z_other))[:sample_size])
                    z_other_sampled = z_other[list(sampler)]
                    
                    # 打印 sample_type_pair[label_mask][high_score_indices[label]] 的值
                    print(f"Label: {label}")
                    print("Values in sample_type_pair[label_mask][high_score_indices[label]]:", sample_type_pair[label_mask][high_score_indices[label]])
    
                    # Compute contrastive losses
                    loss_d1, loss_d2, loss_d3 = compute_losses(z_reference, high_score_features, low_score_features, z_other_sampled)
                    
                    print(f"\tLosses: loss_d1={loss_d1}, loss_d2={loss_d2}, loss_d3={loss_d3}")
                    
                    # 打印 z_reference, high_score_features, low_score_features, z_other_sampled 的形状和内容
#                     print("z_reference:", z_reference)
                    print("\tz_reference shape:", z_reference.shape)

#                     print("high_score_features:", high_score_features)
                    print("\thigh_score_features shape:", high_score_features.shape)

#                     print("low_score_features:", low_score_features)
                    print("\tlow_score_features shape:", low_score_features.shape)

#                     print("z_other_sampled:", z_other_sampled)
                    print("\tz_other_sampled shape:", z_other_sampled.shape)


                    # Add valid losses to total
                    if not torch.isnan(loss_d1) and not torch.isnan(loss_d2) and not torch.isnan(loss_d3):
                        total_contrastive_loss += loss_d1 + loss_d2 + loss_d3
                        valid_contrastive_loss_count += 1
            
            # Average total contrastive loss if there are valid counts
            if valid_contrastive_loss_count > 0:
                total_contrastive_loss /= valid_contrastive_loss_count

            # Compute weighted losses
            weighted_contrastive_loss = total_contrastive_loss * 0.3
            weighted_ucell_loss = total_ucell_loss * 0.3
            weighted_kl_loss = (total_loss_kl1 / self.n_batch) * 10
            weighted_cluster_loss = (total_loss_cluster / self.n_batch) * 1.0
            
            # Compute total loss
            total_loss = weighted_kl_loss + weighted_cluster_loss + weighted_ucell_loss + weighted_contrastive_loss

            print(f'Epoch {epoch + 1}:')
            print(f'  KL Loss: {weighted_kl_loss}')
            print(f'  Cluster Loss: {weighted_cluster_loss}')
            print(f'  UCell Loss: {weighted_ucell_loss}')
            print(f'  Total Contrastive Loss: {weighted_contrastive_loss}')
            print(f'  Total Loss: {total_loss}')
            
            self.gnn_optimizer.zero_grad()
            self.net_optimizer.zero_grad()
            total_loss.backward()
            self.gnn_optimizer.step()
            self.net_optimizer.step()

        return self.gnn

    def train_model(self, indices, RNA_matrix, ini_p1, sample_type, nodes_id):
        self.train()
        print('The training process for the EmitGCL model has started. Please wait.')
        EmitGCL_gnn = self.forward(indices, RNA_matrix, ini_p1, sample_type, nodes_id)
        print('The training for the EmitGCL model has been completed.')
        return EmitGCL_gnn
        
        
def EmitGCL_pred(RNA_matrix, EmitGCL_gnn, indices, nheads, nodes_id, cell_size, device, gene_names, node_dim_reduction_model):
    n_batch = len(indices)
    embedding = []
    l_pre = []
    EmitGCL_result = {}
        
    # Assume global counts of cells and genes are known
    num_cells = RNA_matrix.shape[1]  # Total global cell count
    num_genes = RNA_matrix.shape[0]  # Total global gene count

    # Initialize global attention matrices, creating one matrix for each head
    global_attention_matrices = {f'head-{i+1}': torch.zeros(num_genes, num_cells, device=device) for i in range(nheads)}

    with torch.no_grad():
        for batch_id in tqdm(range(n_batch)):
            gene_index = indices[batch_id]['gene_index']
            cell_index = indices[batch_id]['cell_index']
            gene_feature = torch.tensor(np.array(RNA_matrix[list(gene_index),].todense()), dtype=torch.float32).to(device)
            cell_feature = torch.tensor(np.array(RNA_matrix[:, list(cell_index)].T.todense()), dtype=torch.float32).to(device)
            node_feature = [cell_feature, gene_feature]
            
            gene_cell_sub = RNA_matrix[list(gene_index),][:, list(cell_index)]
            gene_cell_edge_index = torch.LongTensor([np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1], np.nonzero(gene_cell_sub)[1]]).to(device)
            
            edge_index = gene_cell_edge_index
            node_type = torch.LongTensor(np.array(list(np.zeros(len(cell_index))) + list(np.ones(len(gene_index))))).to(device)
            edge_type = torch.LongTensor(np.array(list(np.zeros(gene_cell_edge_index.shape[1])))).to(device)
            
            node_rep = EmitGCL_gnn.forward(node_feature, node_type, edge_index, edge_type).to(device)
            
            # Retrieve current subgraph's attention weights
            attention_weights = node_dim_reduction_model.subgraph_attention_weights[batch_id]
            
            # Get global cell and gene indices for the current subgraph
            global_gene_indices = indices[batch_id]['gene_index']
            global_cell_indices = indices[batch_id]['cell_index']
            
            # Iterate through each edge in the current subgraph and update global attention matrices
            for (local_gene_idx, local_cell_idx), weights_vector in attention_weights.items():
                global_gene_idx = global_gene_indices[local_gene_idx-len(indices[batch_id]['cell_index'])]
                global_cell_idx = global_cell_indices[local_cell_idx]

                # Update global attention matrices for each head
                for head_index in range(nheads):
                    head_key = f'head-{head_index+1}'
                    weight = weights_vector[head_index]
                    # Correctly convert numpy.float32 to torch.FloatTensor and assign
                    global_attention_matrices[head_key][global_gene_idx, global_cell_idx] = torch.tensor(weight, device=device)


            cell_emb = node_rep[node_type == 0]
            gene_emb = node_rep[node_type == 1]

            if device.type == "cuda":
                cell_emb = cell_emb.cpu()
            embedding.append(cell_emb.detach().numpy())

            cell_pre = list(cell_emb.argmax(dim=1).detach().numpy())
            l_pre.extend(cell_pre)

    cell_embedding = np.vstack(embedding)
    cell_clu = np.array(l_pre)
    
    # Create a mapping from the original indices to the shuffled indices
    index_mapping = {original_index: node_id for original_index, node_id in enumerate(list(nodes_id))}

    # Reorder cell_clu
    sorted_cell_clu = [None] * len(cell_clu)
    for original_index, node_id in index_mapping.items():
        sorted_cell_clu[node_id] = cell_clu[original_index]

    # Reorder cell_embedding
    sorted_cell_embedding = np.zeros_like(cell_embedding)
    for original_index, node_id in index_mapping.items():
        sorted_cell_embedding[node_id, :] = cell_embedding[original_index, :]


    EmitGCL_result = {'pred_label': sorted_cell_clu, 'cell_embedding': sorted_cell_embedding, 'attention_weights': global_attention_matrices}
    return EmitGCL_result
    
    
# Save a single attention matrix in .mtx format
def save_attention_matrix_to_mtx(tensor_matrix, save_path, file_name):
    # Convert PyTorch tensor to numpy array
    if isinstance(tensor_matrix, torch.Tensor):
        tensor_matrix = tensor_matrix.cpu().numpy()
    # Convert numpy array to CSR format sparse matrix
    sparse_matrix = csr_matrix(tensor_matrix)
    # Save as .mtx format using mmwrite
    mmwrite(os.path.join(save_path, f"{file_name}.mtx"), sparse_matrix)
    
    
def save_result(EmitGCL_result, output_file, attention_file, save_attention=True):
    """
    Save the results of EmitGCL including attention matrices and other data.
    
    Parameters:
    - EmitGCL_result: dict containing the results of the EmitGCL model, including 'pred_label', 'cell_embedding', and 'attention_weights'.
    - output_file: str, the file path to save the results.
    - attention_file: str, the file path to save the attention matrices.
    - save_attention: bool, whether to save the attention matrices.
    """
    
    # Save attention matrices if save_attention is True
    if save_attention:
        attention_matrices = EmitGCL_result['attention_weights']
        # Iterate through the dictionary and save each attention matrix
        for head, matrix in attention_matrices.items():
            save_attention_matrix_to_mtx(matrix, attention_file, head)
    
    # Save additional results to .npy files
    np.save(output_file + "pred.npy", EmitGCL_result['pred_label'])
    np.save(output_file + "cell_embedding.npy", EmitGCL_result['cell_embedding'])
