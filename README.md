# sarcopenia-network-entropy
**Paper title:** 
**The biphasic and age-dependent impact of Klotho on hallmarks of aging and skeletal muscle function** 

Abstract: 
Aging is accompanied by a disrupted information flow, which results from accumulation of molecular mistakes. These mistakes ultimately give rise to debilitating disorders such as skeletal muscle wasting, or sarcopenia. To estimate the growing “disorderliness” of the aging muscle system, we employed a statistical physics approach to estimate the state parameter, entropy, as a function of genes associated with hallmarks of aging. Although the most prominent structural and functional alterations were observed in the oldest old mice (27-29 months), we found that the escalating network entropy reached an inflection point at old age (22-24 months). To probe the potential for restoration of molecular “order” and reversal of the sarcopenic phenotype, we overexpressed the longevity protein, α-Klotho. Klotho overexpression modulated genes representing all hallmarks of aging in both old and oldest-old mice. However, whereas Klotho improved strength in old mice, intervention failed to induce a benefit beyond the entropic tipping point.

Link to paper: https://www.biorxiv.org/content/10.1101/2020.07.22.207043v1

Link to ShinyApp that classifies genes based on Hallmarks of aging and vice-versa:
https://sruthisivakumar.shinyapps.io/HallmarksAgingGenes/

Link to the network entropy calculator website:
https://network-entropy-calculator.herokuapp.com/

Files and description: 
1. Preprocessing jupyter notebook: creates protein-protein interaction network from RNA seq data
2. Network entropy python script: computes network entropy from nodelist and edgelist
3. Tab files: Raw RNA-seq counts
4. deseq2 and functional enrichment: Differential expression gene lists and functional enrichment results

For further details about mathemtical basis of network entropy, please check GitHub: menicgiulia/Code_network_ensembles
