### Transfer_Learning_Precision_Medicine  
Predictive modeling of drug sensitivity is an important step in precision therapy design and often there is a shortage of suitable data for modeling. Hence we attempt to use data from multiple sources for modeling purposes and the recent advent of large-scale pharmacogenomic studies offers a convenient getaway from this conundrum. Here we use genomic and sensitivity data from two renowned studies - [CCLE](https://portals.broadinstitute.org/ccle) and [GDSC](http://www.cancerrxgene.org/) and combine them through **Transfer Learning (TL)** since samples from different sources cannot be used together directly. We implement two different TL approaches - 
* Latent Variable based Cost Optimization
    * Latent regression prediction
    * Latent-latent prediction
    * Combined latent prediction
* Domain Transfer _via_ Nonlinear Mapping  

The details of these approaches are described in the 2018 paper [Application of transfer learning for cancer drug sensitivity prediction](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2465-y). Below provides an overview of the TL scenarios involved in this implementation. 

![TransferLearningSummary](https://github.com/dhruba018/Transfer_Learning_Precision_Medicine/blob/master/TLsummary.jpg)
