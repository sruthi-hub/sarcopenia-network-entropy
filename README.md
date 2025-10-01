# sarcopenia-network-entropy
**Paper title:** 
**The biphasic and age-dependent impact of Klotho on hallmarks of aging and skeletal muscle function** 

Abstract: 
Aging is accompanied by disrupted information flow, resulting from accumulation of molecular mistakes. These mistakes ultimately give rise to debilitating disorders including skeletal muscle wasting, or sarcopenia. To derive a global metric of growing ‘disorderliness’ of aging muscle, we employed a statistical physics approach to estimate the state parameter, entropy, as a function of genes associated with hallmarks of aging. Escalating network entropy reached an inflection point at old age, while structural and functional alterations progressed into oldest-old age. To probe the potential for restoration of molecular ‘order’ and reversal of the sarcopenic phenotype, we systemically overexpressed the longevity protein, Klotho, via AAV. Klotho overexpression modulated genes representing all hallmarks of aging in old and oldest-old mice, but pathway enrichment revealed directions of changes were, for many genes, age-dependent. Functional improvements were also age-dependent. Klotho improved strength in old mice, but failed to induce benefits beyond the entropic tipping point.

Link to paper: https://elifesciences.org/articles/61138

Link to ShinyApp that classifies genes based on Hallmarks of aging and vice-versa:
https://sruthisivakumar.shinyapps.io/HallmarksAgingGenes/

Link to the network entropy calculator website:
https://network-entropy-calculator.streamlit.app/

Files and description: 
1. Preprocessing jupyter notebook: creates protein-protein interaction network from RNA seq data
2. Network entropy python script: computes network entropy from nodelist and edgelist
3. Tab files: Raw RNA-seq counts
4. deseq2 and functional enrichment: Differential expression gene lists and functional enrichment results

For further details about mathematical basis of network entropy, please check GitHub: menicgiulia/Code_network_ensembles
