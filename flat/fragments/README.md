## Fragments

| Code                                     | Name                          |
|:---------------------------------------- |:------------------------------|
| [CSD](https://www.rcsb.org/ligand/CSD)   | CYSTEINESULFINIC ACID         |
| [DAH](https://www.rcsb.org/ligand/DAH)   | DOPA                          |
| [IAS](https://www.rcsb.org/ligand/IAS)   | ASPARTIC ACID                 |
| [KYN](https://www.rcsb.org/ligand/KYN)   | KYNURENINE                    |
| [OCS](https://www.rcsb.org/ligand/OCS)   | CYSTEINESULFONIC ACID         |
| [OHI](https://www.rcsb.org/ligand/OHI)   |                               |
| [OMT](https://www.rcsb.org/ligand/OMT)   | DIOXYMETHIONINE               |
| [PCA](https://www.rcsb.org/ligand/PCA)   | PYROGLUTAMIC ACID             |
| [SMC](https://www.rcsb.org/ligand/SMC)   | METHYLCYSTEINE                |
| [SME](https://www.rcsb.org/ligand/SME)   | METHIONINE SULFOXIDE          |
| [SNN](https://www.rcsb.org/ligand/SNN)   | SUCCINIMIDE                   |
| [TPO](https://www.rcsb.org/ligand/TPO)   | PHOSPHOTHREONINE              |

## About

```python
from pymol import cmd, chempy

# Save fragment
cmd.load("frag.cif")
model = cmd.get_model(frag)
chempy.io.pkl.toFile(model, "frag.pkl")

# Load fragment
model = chempy.io.pkl.fromFile("frag.pkl")
cmd.load_model(model, "frag")

```