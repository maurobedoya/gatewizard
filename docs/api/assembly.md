# Biological assembly

The coordinates in a PDB or mmCIF file are the asymmetric unit. The biological oligomer or capsid is described by operators already in that file: PDB `REMARK 350` / `BIOMT`, or mmCIF `_pdbx_struct_assembly`. Crystal neighbors from `CRYST1` are not built.

Chain copies are written as mmCIF with names such as `A1` and `A2`. That is the same idea as an official biological-assembly mmCIF (`A`, `A-2`, `A-3`). Classic PDB has a one-letter chain column, so it is not used for the assembly itself. Saving the result as PDB keeps only the first letter of each name.

## Import

```python
from gatewizard.core.assembly import list_assemblies, build_assembly
```

### list_assemblies

```python
list_assemblies(path) -> dict
```

Read the original PDB or mmCIF. Each assembly reports its id, the author or software description, how many chains the expansion produces, and an atom estimate with and without crystal waters.

An assembly whose operators are all the identity is marked `already_biological_unit` and is not expanded. A file with no assembly records returns an empty list and a message that there is nothing to build. Strict NCS (`MTRIX`) is listed only when the file has no biological assembly that moves a chain.

### build_assembly

```python
build_assembly(path, assembly_id, include_waters=False, drop_overlaps=True, dest=None) -> dict
```

Write one assembly as mmCIF. `HOH`, `WAT`, and `DOD` are dropped unless `include_waters` is true. Ligands are kept. Identity operators do not duplicate atoms: the result has `built: false` and says there is nothing to build.

Ions and waters that sit on a symmetry axis are copied once per operator and land on the same coordinates. `drop_overlaps` keeps one hetatm at each of those points. Protein atoms that overlap are reported and left in place. The list entry counts those extras in `stacked_extra_atoms` and `stacked_extra_atoms_with_waters`.

Assemblies of about 80,000 atoms or more set `needs_confirmation` on the list entry so a caller can ask before writing.
