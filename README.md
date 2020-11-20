# Transfer Learning for Precision Medicine  
**Reference:** [Application of transfer learning for cancer drug sensitivity prediction](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2465-y).  

Predictive modeling of drug sensitivity is an important step in precision therapy design and often there is a shortage of suitable data for modeling. Hence we attempt to use data from multiple sources for modeling purposes and the recent advent of large-scale pharmacogenomic studies offers a convenient getaway from this conundrum.  

## Description
### Data
We use *in vitro* transcriptomic and drug sensitivity information from the administration of various anticancer drugs on multiple Cancer cell lines of different subtypes. This data were obtained from two renowned large-scale pharmacogenomic studies - 
   * [Cancer Cell Line Encyclopedia (CCLE)](https://portals.broadinstitute.org/ccle/about/). A collaboration between Broad Institute and Novartis Institutes for Biomedical Research.
   * [Genomics of Drug Sensitivity in Cancer (GDSC)](https://www.cancerrxgene.org/about). A collaboration between Wellcome Sanger Institute (UK) and Center for Molecular Therapeutics at Massachusetts General Hospital Cancer Center (USA). 

There is a considerable overlap between the two studies and we have explored the existing *inconsistencies between datasets* in our previous work in Briefings in Bioinformatics: [Evaluating the consistency of large-scale pharmacogenomic studies](https://academic.oup.com/bib/article-abstract/20/5/1734/5034074). 

### Approaches
We have combined datasets from CCLE and GDSC through **Transfer Learning** since the samples from two different sources cannot be used together directly. To eliminate the distribution shift present in the two sets, we have implemented two different TL approaches - 
#### <ins><b>Latent Variable based Cost Optimization</b></ins>.  
We use the notion of a *Latent variable space* to model the underlying similarities between the genomic and sensitivity datasets and attempt to minimize the discepancies through cost optimization.  
     
![lvco_eqn](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Clarge%20%5Cboldsymbol%7Bw%7D%20%3D%20c_0%20&plus;%20c_p%20%5Cboldsymbol%7Bz%7D_p%20&plus;%20c_s%20%5Cboldsymbol%7Bz%7D_s%20&plus;%20%5Cboldsymbol%7B%5Cvarepsilon%7D%2C%20%5Cqquad%20%5Csum_i%20%7Bc_i%7D%20%3D%201)  
     
where <i><b>z</b><sub>p</sub></i> and <i><b>z</b><sub>s</sub></i> represents the primary (target) and secondary (source) sets, and <i><b>w</b></i> is the underlying latent variable.  
We have implemented three different optimization based approaches - 
  * Latent regression prediction
  * Latent-latent prediction
  * Combined latent prediction
     
#### <ins><b>Domain Transfer _via_ Nonlinear Mapping</b></ins>.  
We have implemented a one-to-one feature mapping between the samples in the primary and secondary datasets using the *Polynomial regression mapping* to transfer the primary data to the secondary space and perform prediction using the larger datasets available in the secondary space. Note that, this approach assumes the existence of a set of *matched samples* between the two sets. 
     
![dtnm_eqn](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Clarge%20%5Cboldsymbol%7Bz%7D_%7Bs%2C%20k%7D%20%3D%20%5Cboldsymbol%7B%5Cgamma%7D_%7Bp%7D%5E%7B%28k%29%7D%20%5Cboldsymbol%7Bz%7D_%7Bp%2C%20k%7D%20&plus;%20%5Cboldsymbol%7B%5Cvarepsilon%7D%5E%7B%28k%29%7D)  
     
where <i><b>z</b><sub>p, k</sub></i> and <i><b>z</b><sub>s, k</sub></i> represents the <i><b>k</b></i>-th feature (gene/drug) in the primary and secondary sets, and <i><b>γ</b></i> is the polynomial mapping coefficient.      

The details of these approaches are described in our 2018 paper: [Application of transfer learning for cancer drug sensitivity prediction](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2465-y). Below provides an overview of the TL scenarios involved in this implementation. 

![TL_summary](https://github.com/dhruba018/Transfer_Learning_Precision_Medicine/blob/master/TLsummary.jpg)

## Cite
If you use either the Domain Transfer TL approach or Latent Variable Cost Optimization TL approach for your research/application, please cite the following paper - 
  > Dhruba, S., Rahman, R. et al. Application of transfer learning for cancer drug sensitivity prediction. BMC Bioinformatics 19, 497 (2018). 
    DOI: https://doi.org/10.1186/s12859-018-2465-y


 you use our work on the exploration of inconsistencies in the large pharmacogenomic studies, please cite the following paper - 
  > Rahman, R., Dhruba, S. R. et al., Evaluating the consistency of large-scale pharmacogenomic studies, Briefings in Bioinformatics, 20 (5), 1734 – 1753 (2019). 
    DOI: https://doi.org/10.1093/bib/bby046
