# EmitGCL Project

## Introduction

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Welcome to the EmitGCL project. This repository contains all the source code for EmitGCL, a graph contrastive learning model designed to detect early metastatic cell groups by integrating knowledge from primary and metastatic sites. EmitGCL excels in identifying subtle cellular differences, with high sensitivity and low false positive rates, especially in rare and early-stage metastatic cells, offering crucial insights for cancer prognosis and treatment.

---

## Source Code Overview

The source code is organized into several files, each handling a key part of the model pipeline:

### 1. `__init__.py`
Marks this directory as a Python package and initializes the module environment.

---

### 2. `conv.py`
- **Purpose**:  
  Defines graph convolutional layers and attention mechanisms for information propagation on cell graphs.  
- **Key Features**:  
  - Supports multiple convolution types, including attention-based variants.  
  - Provides utilities for node feature aggregation and multi-head attention.  

---

### 3. `EmitGCL.py`
- **Purpose**:  
  Core implementation of the **EmitGCL** model.  
- **Key Features**:  
  - Graph construction and node embedding learning.  
  - Contrastive learning framework for identifying early metastatic cells.  
  - Multi-layer graph convolution integration for high-dimensional data.  

---

### 4. `loss_function.py`
- **Purpose**:  
  Implements all loss functions used in training.  
- **Key Features**:  
  - **Contrastive Loss** for positive/negative sample separation.  
  - **PathwayUCell Loss** leveraging pathway-level gene information.  
  - **Label Smoothing** to mitigate class imbalance and overfitting.  

---

### 5. `utils.py`
- **Purpose**:  
  Provides auxiliary functions for data preprocessing and analysis.  
- **Key Features**:  
  - KEGG pathway gene extraction for metastasis-related signals.  
  - Subgraph sampling and clustering initialization.  
  - Data normalization and transformation utilities.  

---

