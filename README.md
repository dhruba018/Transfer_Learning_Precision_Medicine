## Transfer Learning for Precision Medicine  
Predictive modeling of drug sensitivity is an important step in precision therapy design and often there is a shortage of suitable data for modeling. Hence we attempt to use data from multiple sources for modeling purposes and the recent advent of large-scale pharmacogenomic studies offers a convenient getaway from this conundrum. 

### Description
**Data.**  
We use *in vitro* transcriptomic and drug sensitivity information from the administration of various anticancer drugs on multiple Cancer cell lines of different subtypes. This data were obtained from two renowned large-scale pharmacogenomic studies - 
   * [Cancer Cell Line Encyclopedia (CCLE)](https://portals.broadinstitute.org/ccle/about/). A collaboration between Broad Institute and Novartis Institutes for Biomedical Research.
   * [Genomics of Drug Sensitivity in Cancer (GDSC)](https://www.cancerrxgene.org/about). A collaboration between Wellcome Sanger Institute (UK) and Center for Molecular Therapeutics at Massachusetts General Hospital Cancer Center (USA). 

There is a considerable overlap between the two studies and we have explored the existing *inconsistencies between datasets* in our previous work [Evaluating the consistency of large-scale pharmacogenomic studies](https://academic.oup.com/bib/article-abstract/20/5/1734/5034074). 

**Approaches.**  
We have combined datasets from CCLE and GDSC through **Transfer Learning (TL)** since the samples from two different sources cannot be used together directly. To eliminate the distribution shift present in the two sets, we have implemented two different TL approaches - 
   * <ins><b>Latent Variable based Cost Optimization</b></ins>.  
     We use a *Latent variable space* to model the underlying similarities between the genomic and sensitivity datasets and try to minimize the discepancies _via_ cost optimization. We implemented three different subapproaches -      
      * Latent regression prediction
      * Latent-latent prediction
      * Combined latent prediction
   * <ins><b>Domain Transfer _via_ Nonlinear Mapping</b></ins>.  
     We implement a one-to-one sample mapping between primary (target) and secondary (source) datasets using *Polynomial regression mapping*.
     
     ![equation]([img]http://www.sciweavers.org/tex2img.php?eq=%5Cboldsymbol%7Bg%7D_%7Bs%2C%20i%7D%20%3D%20%5Cboldsymbol%7B%5Calpha%7D_p%5E%7B%28i%29%7D%20%5C%2C%20%5Cboldsymbol%7Bg%7D_%7Bp%2C%20i%7D%20%2B%20%5Cboldsymbol%7B%5Cvarepsilon%7D%5E%7B%28i%29%7D%20%5C%5C%0A%5Cboldsymbol%7Bd%7D_%7Bp%2C%20j%7D%20%3D%20%5Cboldsymbol%7B%5Cbeta%7D_s%5E%7B%28j%29%7D%20%5C%2C%20%5Cboldsymbol%7Bd%7D_%7Bs%2C%20j%7D%20%2B%20%5Cboldsymbol%7B%5Cvarepsilon%7D%5E%7B%28j%29%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0[/img])

The details of these approaches are described in the 2018 paper [Application of transfer learning for cancer drug sensitivity prediction](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2465-y). Below provides an overview of the TL scenarios involved in this implementation. 

![TransferLearningSummary](https://github.com/dhruba018/Transfer_Learning_Precision_Medicine/blob/master/TLsummary.jpg)

### Cite
If you use either the Domain Transfer TL approach or Latent Variable Cost Optimization TL approach for your research/application, please cite the following paper - 
> Dhruba, S., Rahman, R., Matlock, K. et al. Application of transfer learning for cancer drug sensitivity prediction. BMC Bioinformatics 19, 497 (2018). 
  DOI: https://doi.org/10.1186/s12859-018-2465-y
