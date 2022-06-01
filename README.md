DSCN: Double target selection guided by CRISPR screening and network model
======================================================================

Introduction
------------

DSCN is a tool that takes multi-omics data and prioritize target combinations that facilitate the development
of novel treatment plan of complex diseases such as cancers.
DSCN derived from the previous work 'SCNrank', which takes the same input and prioritize single targets.
Check <a href="https://link.springer.com/article/10.1186/s12920-020-0681-6">SCNrank</a> if you are interested in our previous work.


Overview
--------
DSCN takes multi-omics data as input from both tissue and cell-line models, construct functional networks for each respectively.
Collect features from tissue network, map them onto cell-line network, where actionable target combinations are predicted.
![Slide1.jpg](https://www.biorxiv.org/content/biorxiv/early/2021/09/06/2021.09.06.459081/F1.large.jpg?width=800&height=600&carousel=1)
**Fig 1.** General workflow of DSCN predicting target combinations.

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
Citation
--------
To be added

License
-------
We follow the MIT License so people are free to check and use the code.

Usage
------
STEP1: Prepare your input

a_file. A tumor expression profile

b_file. A tumor vs normal fold change profile

c_file. A cell-line expression profile

d_file. A cell-line CRISPR screening profile

e_file. A PPI network (STRING human PPI V10 with no cut-offs recommended)

f_file. A target file indicating the designated set of targets.

STEP2: Make a overlapped subset of your input files (a_file'-f_file')

STEP3: Run spectral clustering first

e.g. python DSCN_WG_SC.py a_file' b_file' e_file' > SC_result

STEP4: Score target combinations 

e.g. python DSCN_MP_TN_refine.py f_file' c_file' f_file' e_file' SC_result > ranked_list

STEP5: Estimate subnetwork similarity between tumor and cell-lines (optional)

e.g. python DSCN_sn_similarity.py f_file' a_file' b_file' c_file' d_file' e_file' SC_result 