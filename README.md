## Transfer Learning for Precision Medicine  
Predictive modeling of drug sensitivity is an important step in precision therapy design and often there is a shortage of suitable data for modeling. Hence we attempt to use data from multiple sources for modeling purposes and the recent advent of large-scale pharmacogenomic studies offers a convenient getaway from this conundrum. 

### Description
**Data.**  
We use *in vitro* transcriptomic and drug sensitivity information from the administration of various anticancer drugs on multiple Cancer cell lines of different subtypes. This data were obtained from two renowned large-scale pharmacogenomic studies - 
   * [Cancer Cell Line Encyclopedia (CCLE)](https://portals.broadinstitute.org/ccle/about/). A collaboration between Broad Institute and Novartis Institutes for Biomedical Research.
   * [Genomics of Drug Sensitivity in Cancer (GDSC)](https://www.cancerrxgene.org/about). A collaboration between Wellcome Sanger Institute (UK) and Center for Molecular Therapeutics at Massachusetts General Hospital Cancer Center (USA). 

**Approaches.**  
We combine datasets from CCLE and GDSC through **Transfer Learning (TL)** since samples from different sources cannot be used together directly (we performed detailed comparison for these two studies in our earlier work [in Briefings in Bioinformatics](https://academic.oup.com/bib/article-abstract/20/5/1734/5034074)). We implement two different TL approaches - 
   * Latent Variable based Cost Optimization
      * Latent regression prediction
      * Latent-latent prediction
      * Combined latent prediction
   * Domain Transfer _via_ Nonlinear Mapping  

The details of these approaches are described in the 2018 paper [Application of transfer learning for cancer drug sensitivity prediction](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2465-y). Below provides an overview of the TL scenarios involved in this implementation. 

![TransferLearningSummary](https://github.com/dhruba018/Transfer_Learning_Precision_Medicine/blob/master/TLsummary.jpg)

**Cite.**  
If you use either the Domain Transfer TL approach or Latent Variable Cost Optimization TL approach for your research/application, please cite the following paper - 

      Dhruba, S., Rahman, R., Matlock, K. et al. Application of transfer learning for cancer drug sensitivity prediction. BMC Bioinformatics 19, 497 (2018).               
      https://doi.org/10.1186/s12859-018-2465-y
