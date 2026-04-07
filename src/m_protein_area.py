"""
This Script Computes the projected surface area of a m protein in three orthogonal planes (XY, XZ,
YZ).

We will use this data to estimate the copy number of m proteins per 
virion based on the total surface area of the virion and the packing 
density of m proteins (implemented in protein_analyzer.py).

============================

Author: Syed Mushahid Hussain
"""
from Bio.PDB import MMCIFParser
import numpy as np
from scipy.spatial import ConvexHull

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("4LP7", "4LP7-assembly1.cif")

coords = np.array([a.get_vector().get_array() for a in structure.get_atoms()])

# Project onto each face and compute CONVEX HULL area
# This gives the actual shadow area, not bounding rectangle
faces = [
    ("YZ plane (looking along X)", coords[:, [1, 2]]),
    ("XZ plane (looking along Y)", coords[:, [0, 2]]),
    ("XY plane (looking along Z)", coords[:, [0, 1]])
]

print("ACTUAL PROJECTED AREAS (Convex Hull)")
print("=" * 50)
for name, proj in faces:
    hull = ConvexHull(proj)
    area_A2 = hull.volume  # In 2D, ConvexHull.volume = area
    area_nm2 = area_A2 / 100  # Convert Å² to nm²
    print(f"{name}:")
    print(f"  Convex hull area: {area_nm2:.1f} nm²")
    print()

# Compare with bounding box approach
d = sorted(coords.max(axis=0) - coords.min(axis=0), reverse=True)
bbox_face = d[0] * d[1] / 100
print(f"Bounding box largest face: {bbox_face:.1f} nm²")
print(f"Bounding box × 0.65:       {bbox_face * 0.65:.1f} nm²")
print()
print("The convex hull gives the TRUE projected area")
print("without needing a packing fraction assumption")