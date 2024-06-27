# Comparative Systems Service

## Overview

The Comparative Systems Service combines together functionality from three different tools and data from the legacy PATRIC system: the **Protein Family Sorter**, **Pathway Comparison Tool**, and **Subsystems Data**. Subsystems (1,2), are functional roles that together implement a specific biological process or structural complex and can also be generalized as pathways. Up to 500 genomes can be compared. All three tools use the two protein families, PATtyFams (3), that are assigned in the BV-BRC annotation process known as RASTtk (4). The global families, known as PGFams, can be used for cross genus comparisons. The local families, PLFams, are for intra-genus comparisons. Pathway maps are represented using KEGG (5).



## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

This module provides the following application specfication(s):
* [ComparativeSystems](app_specs/ComparativeSystems.md)


## See also

* [Comparative Systems Service Quick Reference](https://www.bv-brc.org/docs/quick_references/services/comparative_systems.html)
  * [Comparative Systems Service](https://www.bv-brc.org/docs/https://www.bv-brc.org/app/ComparativeSystems.html)
  * [Comparative Systems Service Tutorial](https://www.bv-brc.org/docs//tutorial/comparative_systems/comparative_systems.html)
  * [Pathway Comparison Tool Quick Reference Guide](https://www.bv-brc.org/docs//quick_references/other/pathway_comparison_tool.html)
  * [Protein Family Sorter Quick Reference Guide](https://www.bv-brc.org/docs//quick_references/other/protein_family_sorter.html)
  * [Subystems Data Quick Reference Guide](https://www.bv-brc.org/docs//quick_references/other/subsystems_data.html) 



## References

1. Overbeek, R. et al. The subsystems approach to genome annotation and its use in the project to annotate 1000 genomes. Nucleic acids research 33, 5691-5702 (2005).
2. Overbeek, R. et al. The SEED and the Rapid Annotation of microbial genomes using Subsystems Technology (RAST). 42, D206-D214 (2013).
3. Davis, J. J. et al. PATtyFams: Protein families for the microbial genomes in the PATRIC database. 7, 118 (2016).
4. Brettin, T. et al. RASTtk: a modular and extensible implementation of the RAST algorithm for building custom annotation pipelines and annotating batches of genomes. Scientific reports 5, 8365 (2015).
5. Kanehisa, M., Furumichi, M., Sato, Y., Kawashima, M. & Ishiguro-Watanabe, M. KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Research (2022).

