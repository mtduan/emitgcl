import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from scipy.stats import zscore
from collections import defaultdict


# Define the similarity function
def sim(z_i, z_j):
    # Normalize the feature vectors and calculate cosine similarity
    z_i_norm = F.normalize(z_i, p=2, dim=1)
    z_j_norm = F.normalize(z_j, p=2, dim=1)
    return torch.mm(z_i_norm, z_j_norm.t())

# Define the contrastive loss function
def contrastive_loss(z_i, z_j, z_k, temperature):
    # Calculate similarities for the positive and negative pairs
    pos_sim = sim(z_i, z_j) / temperature
    neg_sim = sim(z_i, z_k) / temperature

    # Compute the loss
    exp_pos_sim = torch.exp(pos_sim)
    exp_neg_sim = torch.exp(neg_sim).sum(dim=1, keepdim=True)  # Sum over the negative samples
    loss = -torch.log(exp_pos_sim / (exp_pos_sim + exp_neg_sim))

    return loss.mean()

def compute_losses(z_lymph, high_metastatic_features, low_metastatic_features, z_other_sampled,
                   temperature_d1=0.1, temperature_d2=0.5, temperature_d3=1.0):
    loss_d1 = contrastive_loss(z_lymph, z_lymph, high_metastatic_features, temperature_d1)
    loss_d2 = contrastive_loss(high_metastatic_features, z_lymph, low_metastatic_features, temperature_d2)
    loss_d3 = contrastive_loss(low_metastatic_features, high_metastatic_features, z_other_sampled, temperature_d3)

    return loss_d1, loss_d2, loss_d3

class PathwayUCellLoss(nn.Module):
    def __init__(self, pathway_genes, gene_name_to_index, device):
        super(PathwayUCellLoss, self).__init__()
        self.pathway_genes = pathway_genes
        self.gene_name_to_index = gene_name_to_index
        self.device = device

    def forward(self, RNA_matrix, labels, sample_type):
        pathway_indices = {pathway: [self.gene_name_to_index[gene] for gene in genes 
                                     if gene in self.gene_name_to_index] 
                           for pathway, genes in self.pathway_genes.items()}

        pathway_scores = defaultdict(list)
        for pathway, indices in pathway_indices.items():
            if indices:
                pathway_expression = RNA_matrix[:, indices]
                pathway_expression_tensor = torch.tensor(pathway_expression.todense(), dtype=torch.float32).to(self.device)
                pathway_scores[pathway] = pathway_expression_tensor.mean(dim=1)

        unique_labels, counts = np.unique(labels.cpu(), return_counts=True)
        underrepresented_labels = unique_labels[counts < 0.05 * len(labels)]
        ucell_loss = 0
        high_metastatic_indices = {}
        low_metastatic_indices = {}
        for label in underrepresented_labels:
            label_mask = np.array((labels == label).cpu())
            tumor_mask = (np.array(sample_type) == 'Tumor') & label_mask
            tumor_mask_label = tumor_mask[label_mask]
            total_pathway_scores = torch.sum(torch.stack([pathway_scores[pathway][label_mask] 
                                                          for pathway in pathway_scores]), dim=0)
            z_scores = zscore(total_pathway_scores.cpu().numpy())
            high_metastatic_indices[label] = np.where((z_scores > 1.64) & tumor_mask_label)[0]
            low_metastatic_indices[label] = np.where((z_scores <= 1.64) & tumor_mask_label)[0]
            ucell_loss += torch.mean(total_pathway_scores)
        
        ucell_loss_tensor = -torch.tensor(ucell_loss, dtype=torch.float32) / len(underrepresented_labels) if underrepresented_labels.size > 0 else torch.tensor(0.0, dtype=torch.float32)
        return ucell_loss_tensor, high_metastatic_indices, low_metastatic_indices
        
        
class LabelSmoothing(nn.Module):
    """NLL loss with label smoothing.
    """

    def __init__(self, smoothing=0.0):
        """Constructor for LabelSmoothing module.
        :param smoothing: Label smoothing factor
        """
        super(LabelSmoothing, self).__init__()
        self.confidence = 1.0 - smoothing
        self.smoothing = smoothing

    def forward(self, x, target):
        logprobs = torch.nn.functional.log_softmax(x, dim=-1)
        nll_loss = -logprobs.gather(dim=-1, index=target.unsqueeze(1))
        nll_loss = nll_loss.squeeze(1)
        smooth_loss = -logprobs.mean(dim=-1)
        loss = self.confidence * nll_loss + self.smoothing * smooth_loss
        return loss.mean()