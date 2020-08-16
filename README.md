### Transfer_Learning_Precision_Medicine  
Predictive modeling of drug sensitivity is an important step in precision therapy design and often there is a shortage of suitable data for modeling. Here we use data from two publicly available large pharmacogenomic studies - [CCLE](https://portals.broadinstitute.org/ccle) and [GDSC](http://www.cancerrxgene.org/) for modeling. Since data from different sources cannot be used together directly, we implement two transfer learning approaches -  
* Latent Variable Cost Optimization
* Domain Transfer via Nonlinear Mapping  
Details of the approaches are described in [this 2018 paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2465-y). 
