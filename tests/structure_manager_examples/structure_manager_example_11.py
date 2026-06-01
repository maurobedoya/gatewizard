import os
import tempfile
from gatewizard.core.structure_manager import StructureManager, StructureError

viewer = StructureManager()

# Backbone atoms (N, CA, C, O) from residues 1-20 of 2MVJ.
# Contains an alpha-helix (residues 5-19) that psique can detect.
pdb_content = """\
HEADER    BACKBONE OF 2MVJ RESIDUES 1-20
ATOM      1  N   MET A   1      12.463  -9.010  13.613  1.00  1.00           N
ATOM      2  CA  MET A   1      12.443  -7.552  13.317  1.00  1.00           C
ATOM      3  C   MET A   1      11.302  -6.897  14.099  1.00  1.00           C
ATOM      4  O   MET A   1      10.705  -7.511  14.982  1.00  1.00           O
ATOM      9  N   LYS A   2      11.011  -5.641  13.765  1.00  1.00           N
ATOM     10  CA  LYS A   2       9.947  -4.887  14.428  1.00  1.00           C
ATOM     11  C   LYS A   2       8.585  -5.480  14.086  1.00  1.00           C
ATOM     12  O   LYS A   2       7.578  -5.120  14.696  1.00  1.00           O
ATOM     18  N   PHE A   3       8.562  -6.358  13.080  1.00  1.00           N
ATOM     19  CA  PHE A   3       7.325  -6.995  12.609  1.00  1.00           C
ATOM     20  C   PHE A   3       6.359  -5.955  12.027  1.00  1.00           C
ATOM     21  O   PHE A   3       5.390  -6.304  11.351  1.00  1.00           O
ATOM     29  N   TYR A   4       6.650  -4.674  12.270  1.00  1.00           N
ATOM     30  CA  TYR A   4       5.852  -3.563  11.761  1.00  1.00           C
ATOM     31  C   TYR A   4       6.211  -3.358  10.285  1.00  1.00           C
ATOM     32  O   TYR A   4       5.510  -2.668   9.545  1.00  1.00           O
ATOM     41  N   THR A   5       7.302  -3.996   9.874  1.00  1.00           N
ATOM     42  CA  THR A   5       7.765  -3.930   8.497  1.00  1.00           C
ATOM     43  C   THR A   5       6.719  -4.546   7.560  1.00  1.00           C
ATOM     44  O   THR A   5       6.375  -3.963   6.532  1.00  1.00           O
ATOM     48  N   ILE A   6       6.190  -5.718   7.941  1.00  1.00           N
ATOM     49  CA  ILE A   6       5.155  -6.387   7.142  1.00  1.00           C
ATOM     50  C   ILE A   6       3.900  -5.510   7.126  1.00  1.00           C
ATOM     51  O   ILE A   6       3.257  -5.341   6.090  1.00  1.00           O
ATOM     56  N   LYS A   7       3.556  -4.976   8.296  1.00  1.00           N
ATOM     57  CA  LYS A   7       2.368  -4.139   8.440  1.00  1.00           C
ATOM     58  C   LYS A   7       2.417  -2.939   7.502  1.00  1.00           C
ATOM     59  O   LYS A   7       1.446  -2.649   6.802  1.00  1.00           O
ATOM     65  N   LEU A   8       3.537  -2.227   7.513  1.00  1.00           N
ATOM     66  CA  LEU A   8       3.683  -1.040   6.683  1.00  1.00           C
ATOM     67  C   LEU A   8       3.602  -1.416   5.212  1.00  1.00           C
ATOM     68  O   LEU A   8       2.948  -0.751   4.419  1.00  1.00           O
ATOM     73  N   ALA A   9       4.262  -2.505   4.854  1.00  1.00           N
ATOM     74  CA  ALA A   9       4.249  -2.968   3.470  1.00  1.00           C
ATOM     75  C   ALA A   9       2.803  -3.133   2.979  1.00  1.00           C
ATOM     76  O   ALA A   9       2.452  -2.708   1.880  1.00  1.00           O
ATOM     78  N   LYS A  10       1.968  -3.727   3.829  1.00  1.00           N
ATOM     79  CA  LYS A  10       0.551  -3.916   3.507  1.00  1.00           C
ATOM     80  C   LYS A  10      -0.150  -2.562   3.402  1.00  1.00           C
ATOM     81  O   LYS A  10      -0.936  -2.325   2.485  1.00  1.00           O
ATOM     87  N   PHE A  11       0.147  -1.681   4.350  1.00  1.00           N
ATOM     88  CA  PHE A  11      -0.448  -0.340   4.375  1.00  1.00           C
ATOM     89  C   PHE A  11      -0.163   0.373   3.058  1.00  1.00           C
ATOM     90  O   PHE A  11      -1.068   0.881   2.404  1.00  1.00           O
ATOM     98  N   LEU A  12       1.106   0.371   2.669  1.00  1.00           N
ATOM     99  CA  LEU A  12       1.514   0.988   1.405  1.00  1.00           C
ATOM    100  C   LEU A  12       0.887   0.229   0.237  1.00  1.00           C
ATOM    101  O   LEU A  12       0.425   0.834  -0.728  1.00  1.00           O
ATOM    106  N   GLY A  13       0.891  -1.100   0.324  1.00  1.00           N
ATOM    107  CA  GLY A  13       0.333  -1.932  -0.737  1.00  1.00           C
ATOM    108  C   GLY A  13      -1.169  -1.735  -0.857  1.00  1.00           C
ATOM    109  O   GLY A  13      -1.790  -2.148  -1.836  1.00  1.00           O
ATOM    110  N   GLY A  14      -1.757  -1.091   0.145  1.00  1.00           N
ATOM    111  CA  GLY A  14      -3.194  -0.806   0.137  1.00  1.00           C
ATOM    112  C   GLY A  14      -3.420   0.524  -0.569  1.00  1.00           C
ATOM    113  O   GLY A  14      -4.313   0.667  -1.404  1.00  1.00           O
ATOM    114  N   ILE A  15      -2.571   1.493  -0.221  1.00  1.00           N
ATOM    115  CA  ILE A  15      -2.618   2.835  -0.802  1.00  1.00           C
ATOM    116  C   ILE A  15      -2.152   2.813  -2.253  1.00  1.00           C
ATOM    117  O   ILE A  15      -2.687   3.540  -3.088  1.00  1.00           O
ATOM    122  N   VAL A  16      -1.114   2.020  -2.551  1.00  1.00           N
ATOM    123  CA  VAL A  16      -0.570   1.990  -3.912  1.00  1.00           C
ATOM    124  C   VAL A  16      -1.698   1.740  -4.925  1.00  1.00           C
ATOM    125  O   VAL A  16      -1.785   2.408  -5.956  1.00  1.00           O
ATOM    129  N   ARG A  17      -2.576   0.798  -4.587  1.00  1.00           N
ATOM    130  CA  ARG A  17      -3.730   0.477  -5.416  1.00  1.00           C
ATOM    131  C   ARG A  17      -4.695   1.668  -5.437  1.00  1.00           C
ATOM    132  O   ARG A  17      -5.274   1.987  -6.472  1.00  1.00           O
ATOM    140  N   ALA A  18      -4.880   2.304  -4.276  1.00  1.00           N
ATOM    141  CA  ALA A  18      -5.800   3.442  -4.169  1.00  1.00           C
ATOM    142  C   ALA A  18      -5.396   4.577  -5.114  1.00  1.00           C
ATOM    143  O   ALA A  18      -6.249   5.205  -5.738  1.00  1.00           O
ATOM    145  N   MET A  19      -4.092   4.828  -5.217  1.00  1.00           N
ATOM    146  CA  MET A  19      -3.578   5.878  -6.081  1.00  1.00           C
ATOM    147  C   MET A  19      -3.933   5.591  -7.539  1.00  1.00           C
ATOM    148  O   MET A  19      -4.546   6.409  -8.217  1.00  1.00           O
ATOM    153  N   LEU A  20      -3.542   4.417  -8.002  1.00  1.00           N
ATOM    154  CA  LEU A  20      -3.823   3.999  -9.374  1.00  1.00           C
ATOM    155  C   LEU A  20      -5.327   3.888  -9.559  1.00  1.00           C
ATOM    156  O   LEU A  20      -5.879   4.265 -10.593  1.00  1.00           O
END
"""

with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)

    # Default SS (assigned automatically at load time)
    print("SS after load (auto):")
    print(f"  {viewer.get_secondary_structure_summary()}")

    # Reassign using the heuristic method
    ss = viewer.assign_secondary_structure("heuristic")
    print(f"SS after heuristic: {ss}")

    # Reassign from PDB HELIX/SHEET records (none in this file)
    try:
        ss = viewer.assign_secondary_structure("pdb_records")
        print(f"SS after pdb_records: {ss}")
    except StructureError as e:
        print(f"pdb_records: {e}")

    # Try psique (may not be installed)
    try:
        ss = viewer.assign_secondary_structure("psique")
        print(f"SS after psique: {ss}")
    except StructureError as e:
        print(f"psique not available: {e}")

    # Auto method (same priority as load)
    ss = viewer.assign_secondary_structure("auto")
    print(f"SS after auto: {ss}")
finally:
    os.unlink(tmp_path)
