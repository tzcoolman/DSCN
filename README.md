DSCN: Double target selection guided by CRISPR screening and network model
======================================================================

Introduction
------------

DSCN is a tool that takes multi-omics data and prioritize target combinations that facilitate the development
of novel treatment plan of complex diseases such as cancers.
DSCN is derived from the previous work 'SCNrank', which takes the same input and prioritize single targets.
Check <a href="https://link.springer.com/article/10.1186/s12920-020-0681-6">SCNrank</a> if you are interested in our previous work.


Overview
--------
DSCN takes multi-omics data as input from both tissue and cell-line models, construct functional networks for each respectively.
Collect features from tissue network, map them onto cell-line network, where actionable target combinations are predicted.
![Slide1.jpg](https://www.biorxiv.org/content/biorxiv/early/2021/09/06/2021.09.06.459081/F1.large.jpg?width=800&height=600&carousel=1)
Liu E, Wu X, Wang L, et al. 2021.

**Fig 1.** General workflow of DSCN predicting target combinations.
DSCN generally contains the following steps:
1. Building Tissue- and Cell-line- specific functional networks.
2. Obtain molecular features from the tissue network via spectral clustering.
3. Map features onto the cell-line network.
4. Select first target from the cell-line network.
5. Select second target (target combination) based on the adjusted cell-line network given the first target.
6. (Optional) add P-values on top the selected target combination.

Quickstart
----------
Run clustering on your tissue data as:



Installation
----------
DSCN has been developed under python/2.7.16 but supports python3.

DSCN doesnot require additional installation. It has the minimum requirements stated as:

numpy>=1.9
multipleprocessing
sklean>=0.16
networkx>=2.0

These packages can be install via pip:
```
pip install numpy==1.6
pip install networkx==2.0
...
```
Or via bioconda:
```
conda install -c bioconda networkx
...
```
Or manually by downloading the packages, installing them and setting up environment variables accordingly:
```
python setup.py install
...
```
Citation
--------
To be added


Input files
-----------
To run DSCN, input files need to be prepared, including:

a_file. A tumor-tissue **expression profile**. It can be sequencing data or array expression data.

b_file. A tumor vs normal tissue **fold change (FC)** profile
**Fold change is a value between 0 and positive infinity. It can not be log-fold-change (logFC). logFC would have negative value that would generate negative [eigen-value](https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors) and result in failure of the spectral clustering[1].**

c_file. A cell-line **expression profile**. It can be sequencing data or array expression data.

d_file. A cell-line **CRISPR screening** profile. CRISPR screening data of pancreatic cell-lines from [DepMap](https://depmap.org/portal/) is used in the example but RNAi data is also supported.

e_file. A Protein-protein interaction (**PPI**) network. [STRING](https://string-db.org/cgi/download?sessionId=bZvjml3QVCPG) human protein links V10 with no cut-off is used as the example. PPI networks from other sources such as PathPPI can also be used.

f_file. A target file indicating the designated set of targets. A subset of target genes of FDA approved drug from [DrugBank](https://go.drugbank.com/) database is used in the example. But any gene (target) can be potentially used.

**Genesets in a_file, b_file, c_file, d_file and e_file need to be identical. Genes in f_file need to be a subset of the identical geneset.**

A script in directory toolbox can be used for creating a consistent subset of a-f files.
```
python find_overlap.py a_file b_file c_file d_file f_file e_file output_path
```
a subset of a-f files will be created under the output_path. 

Usage
------
STEP1: Prepare your input

STEP2: Run spectral clustering first

e.g. python DSCN_WG_SC.py a_file' b_file' e_file' > SC_result

STEP3: Score target combinations 

e.g. python DSCN_MP_TN_refine.py f_file' c_file' f_file' e_file' SC_result > ranked_list

STEP4: Estimate subnetwork similarity between tumor and cell-lines (optional)

e.g. python DSCN_sn_similarity.py f_file' a_file' b_file' c_file' d_file' e_file' SC_result

References
----------
[1] Von Luxburg U. A tutorial on spectral clustering[J]. Statistics and computing, 2007, 17(4): 395-416.
[2] 

License
-------
DSCN follows the MIT License.