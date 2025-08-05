```{eval-rst}
.. currentmodule:: pertpy
```

# Metadata

The metadata module provides tooling to annotate perturbations by querying databases.
Such metadata can aid with the development of biologically informed models and can be used for enrichment tests.

## Cell line

This module allows for the retrieval of various types of information related to cell lines,
including cell line annotation, bulk RNA and protein expression data.

Available databases for cell line metadata:

-   [The Cancer Dependency Map Project at Broad](https://depmap.org/portal/)
-   [The Cancer Dependency Map Project at Sanger](https://depmap.sanger.ac.uk/)
-   [Genomics of Drug Sensitivity in Cancer (GDSC)](https://www.cancerrxgene.org/)

## Compound

The Compound module enables the retrieval of various types of information related to compounds of interest, including the most common synonym, pubchemID and canonical SMILES.

Available databases for compound metadata:

-   [PubChem](https://pubchem.ncbi.nlm.nih.gov/)

## Mechanism of Action

This module aims to retrieve metadata of mechanism of action studies related to perturbagens of interest, depending on the molecular targets.

Available databases for mechanism of action metadata:

-   [CLUE](https://clue.io/)

## Drug

This module allows for the retrieval of Drug target information.

Available databases for drug metadata:

-   [chembl](https://www.ebi.ac.uk/chembl/)

```{eval-rst}
.. autosummary::
    :toctree: metadata
    :recursive:

    metadata.CellLine
    metadata.Compound
    metadata.Moa
    metadata.Drug
    metadata.LookUp
```
