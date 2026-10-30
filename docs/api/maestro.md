# Maestro properties

Read connection-table (CT) properties from Schrödinger Maestro `.mae` / `.maegz` files: Glide docking scores, MM/GBSA components, IFD / Prime energies, and other `r_` / `i_` / `b_` / `s_` fields.

The desktop Visualize page uses the same functions for the entry-picker score column, the headline under each loaded pose, and the Maestro properties table. Energies are reported in **kcal/mol**, as in Schrödinger.

Requires RDKit with `MaeMolSupplier`. Open Babel alone cannot supply these scores.

## Import

```python
from gatewizard.tools.maestro import (
    list_maestro_entries,
    read_maestro_properties,
    pick_headline_property,
    is_maestro_path,
    MaestroError,
)
```

### list_maestro_entries

```python
list_maestro_entries(path) -> dict
```

Walk every connection table in a Maestro file. Returns:

| Key | Meaning |
| --- | --- |
| `sourcePath` | Absolute path of the file that was read |
| `entries` | One dict per CT (see below) |

Each entry:

| Key | Meaning |
| --- | --- |
| `index` | Supplier index (0-based). Use this when selecting CTs. |
| `title` | CT title (`_Name` / `s_m_title`), or `None` |
| `label` | `title · filename` or `CT N · filename` |
| `atomCount` | Number of atoms in the CT |
| `kind` | Always `"maestro_ct"` |
| `properties` | List of `{ "name", "type", "value" }` |

Property `type` is `real`, `integer`, `boolean`, or `string`. Private RDKit names (leading `_`) are omitted. String values longer than 500 characters are omitted (Maestro sometimes stores an embedded PDB in a property).

Raises `FileNotFoundError` when the path is missing, and `MaestroError` when the path is not a Maestro file or RDKit cannot parse it.

### read_maestro_properties

```python
read_maestro_properties(path, indices=None) -> list[dict]
```

Same CT property bags as `list_maestro_entries`, optionally limited to supplier `indices`. Omitting `indices` returns every CT. Each list item has `index`, `title`, `label`, `atomCount`, and `properties`.

Raises `ValueError` when an index is not in the file.

### pick_headline_property

```python
pick_headline_property(properties) -> dict | None
```

Chooses one real score for a short summary, in this order:

1. `r_psp_MMGBSA_dG_Bind` (MMGBSA ΔG bind)
2. `r_i_docking_score` (docking score)
3. `r_i_glide_gscore` (Glide gscore)
4. `r_psp_Prime_Energy` (Prime energy)
5. the first other real property

Returns `None` when there is no real property. The Visualize entry picker and structure rows use this order for **Auto**. In the app you can pick another property with **Show score** (picker) or **On structure rows** (properties table); that choice is stored in the browser and is not part of this Python helper.

### is_maestro_path

```python
is_maestro_path(path) -> bool
```

True for `.mae`, `.maegz`, and `.mae.gz`.

### maestro_ct_properties

```python
maestro_ct_properties(mol) -> list[dict]
```

Extract CT properties from an RDKit molecule that already came from a Maestro supplier. Prefer `list_maestro_entries` or `read_maestro_properties` when you only have a file path.

## Example

```python
from gatewizard.tools.maestro import (
    list_maestro_entries,
    pick_headline_property,
    read_maestro_properties,
)

info = list_maestro_entries("poses.maegz")
for entry in info["entries"]:
    head = pick_headline_property(entry["properties"])
    if head:
        print(entry["label"], head["name"], head["value"])

# Only the first two CTs
rows = read_maestro_properties("poses.maegz", indices=[0, 1])
bind = next(
    (p["value"] for p in rows[0]["properties"] if p["name"] == "r_psp_MMGBSA_dG_Bind"),
    None,
)
```

## Desktop app

On Visualize, **Open** a Maestro file, pick CTs, then use the list-ordered toolbar button for the full property table. Scores are kept in memory for the session. Opening a saved viewpoint or animation re-reads them from the original `.mae` / `.maegz`; they are not stored inside those project files.
