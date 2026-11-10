"""Build LAMMPS data files with red cells, DPD solvent and frozen walls.

Atom types used by all examples:
  1  membrane vertices (each cell is one molecule)
  2  solvent particles
  3  wall particles (frozen)

Only numpy is needed.  The examples import this file through
sys.path, so it can stay in tools/.

If you use this code, please cite Z. Chai et al., Biophys. J. (2026),
https://doi.org/10.1016/j.bpj.2026.10.021.
"""

import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MESH = os.path.join(HERE, "..", "mesh", "rbc500.txt")
INF = 1e20    # unbounded block face


def read_mesh(path=MESH):
    """Return vertices (N, 3) and triangles (M, 3, 0-based) of a mesh file."""
    lines = [l.split("#")[0].strip() for l in open(path)]
    lines = [l for l in lines if l]
    nv = int(lines[0].split()[0])
    nt = int(lines[1].split()[0])
    iv = lines.index("Vertices") + 1
    it = lines.index("Triangles") + 1
    verts = np.array([[float(v) for v in l.split()[1:4]] for l in lines[iv:iv + nv]])
    tris = np.array([[int(v) - 1 for v in l.split()[1:4]] for l in lines[it:it + nt]])
    return verts, tris


def topology(tris):
    """Edges and edge-sharing triangle pairs of a closed, oriented mesh.

    A dihedral (i1, i2, i3, i4) holds the triangles i1-i2-i3 and i4-i3-i2,
    both counterclockwise seen from outside, as dihedral style bend expects.
    """
    third = {}
    for a, b, c in tris:
        for p, q, r in ((a, b, c), (b, c, a), (c, a, b)):
            if (p, q) in third:
                raise ValueError("edge %d-%d is used twice in the same direction" % (p, q))
            third[(p, q)] = r
    bonds, dihedrals = [], []
    for (p, q), r in third.items():
        if p > q:
            continue
        if (q, p) not in third:
            raise ValueError("edge %d-%d has only one triangle: the mesh is not closed" % (p, q))
        bonds.append((p, q))
        dihedrals.append((r, p, q, third[(q, p)]))
    return np.array(bonds), np.array(dihedrals)


def area_volume(x, tris):
    """Total area and enclosed volume (positive for outward triangles)."""
    p = x[tris]
    n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    area = 0.5 * np.linalg.norm(n, axis=1).sum()
    volume = np.einsum("ij,ij->i", n, p.sum(axis=1)).sum() / 18.0
    return area, volume


def rotation(axis, degrees):
    """Rotation matrix about the x, y or z axis."""
    t = np.radians(degrees)
    c, s = np.cos(t), np.sin(t)
    i = "xyz".index(axis)
    j, k = (i + 1) % 3, (i + 2) % 3
    m = np.eye(3)
    m[j, j], m[j, k], m[k, j], m[k, k] = c, -s, s, c
    return m


def inside(points, blocks):
    """True for points inside any block (xlo, xhi, ylo, yhi, zlo, zhi)."""
    hit = np.zeros(len(points), dtype=bool)
    for b in blocks:
        lo, hi = np.array(b[0::2]), np.array(b[1::2])
        hit |= np.all((points > lo) & (points < hi), axis=1)
    return hit


class System:
    """Atoms and membrane topology of one simulation box."""

    def __init__(self, box, seed=1):
        self.box = np.array(box, dtype=float).reshape(3, 2)
        self.rng = np.random.default_rng(seed)
        self.x = []
        self.type = []
        self.mol = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.ncell = 0
        self.natoms = 0

    def add_cell(self, center, rot=None, mesh=MESH):
        """Add one cell; rot is a 3x3 matrix applied to the mesh before shifting."""
        verts, tris = read_mesh(mesh)
        bonds, dihedrals = topology(tris)
        if rot is not None:
            verts = verts @ np.asarray(rot).T
        verts = verts - verts.mean(axis=0) + np.asarray(center, dtype=float)
        lo, hi = self.box[:, 0], self.box[:, 1]
        if np.any(verts <= lo) or np.any(verts >= hi):
            raise ValueError("cell at %s does not fit in the box" % (center,))
        self.ncell += 1
        off = self.natoms + 1
        self.x.append(verts)
        self.type.append(np.full(len(verts), 1))
        self.mol.append(np.full(len(verts), self.ncell))
        self.bonds.append(bonds + off)
        self.angles.append(tris + off)
        self.dihedrals.append(dihedrals + off)
        self.natoms += len(verts)
        return verts

    def add_particles(self, atype, density, region=None, exclude=()):
        """Fill a block with randomly placed particles at the given density.

        region is (xlo, xhi, ylo, yhi, zlo, zhi), clipped to the box; points
        inside any block of exclude are left out.
        """
        r = self.box.copy()
        if region is not None:
            r[:, 0] = np.maximum(r[:, 0], np.array(region[0::2], dtype=float))
            r[:, 1] = np.minimum(r[:, 1], np.array(region[1::2], dtype=float))
        if np.any(r[:, 1] <= r[:, 0]):
            raise ValueError("empty region %s" % (region,))
        # the excluded volume is estimated from a large sample, then the
        # region is filled with that many points that miss the excluded blocks
        probe = self.rng.uniform(r[:, 0], r[:, 1], size=(200000, 3))
        frac = 1.0 - inside(probe, exclude).mean() if exclude else 1.0
        n = int(round(density * np.prod(r[:, 1] - r[:, 0]) * frac))
        pts = np.empty((0, 3))
        while len(pts) < n:
            p = self.rng.uniform(r[:, 0], r[:, 1], size=(2 * n + 100, 3))
            if exclude:
                p = p[~inside(p, exclude)]
            pts = np.vstack([pts, p])
        pts = pts[:n]
        self.x.append(pts)
        self.type.append(np.full(n, atype))
        self.mol.append(np.zeros(n, dtype=int))
        self.natoms += n
        return n

    def write(self, path, title):
        x = np.vstack(self.x)
        atype = np.concatenate(self.type)
        mol = np.concatenate(self.mol)
        bonds = np.vstack(self.bonds) if self.bonds else np.empty((0, 2), int)
        angles = np.vstack(self.angles) if self.angles else np.empty((0, 3), int)
        dihedrals = np.vstack(self.dihedrals) if self.dihedrals else np.empty((0, 4), int)
        with open(path, "w") as f:
            f.write("%s\n\n" % title)
            f.write("%d atoms\n3 atom types\n" % len(x))
            f.write("%d bonds\n1 bond types\n" % len(bonds))
            f.write("%d angles\n1 angle types\n" % len(angles))
            f.write("%d dihedrals\n1 dihedral types\n\n" % len(dihedrals))
            for k, ax in enumerate("xyz"):
                f.write("%.6f %.6f %slo %shi\n" % (self.box[k, 0], self.box[k, 1], ax, ax))
            f.write("\nMasses\n\n1 1.0\n2 1.0\n3 1.0\n")
            f.write("\nAtoms # molecular\n\n")
            for i in range(len(x)):
                f.write("%d %d %d %.6f %.6f %.6f\n" % (i + 1, mol[i], atype[i], *x[i]))
            if len(bonds):
                f.write("\nBonds\n\n")
                for i, b in enumerate(bonds):
                    f.write("%d 1 %d %d\n" % (i + 1, *b))
                f.write("\nAngles\n\n")
                for i, a in enumerate(angles):
                    f.write("%d 1 %d %d %d\n" % (i + 1, *a))
                f.write("\nDihedrals\n\n")
                for i, d in enumerate(dihedrals):
                    f.write("%d 1 %d %d %d %d\n" % (i + 1, *d))
        counts = {t: int((atype == t).sum()) for t in (1, 2, 3)}
        print("%s: %d atoms (membrane %d, solvent %d, wall %d), %d cells"
              % (path, len(x), counts[1], counts[2], counts[3], self.ncell))
        return counts


def check_mesh(path=MESH, area=132.8686, volume=92.4523, tol=1e-3):
    """Stop if the mesh is not the closed, outward-oriented cell we expect."""
    verts, tris = read_mesh(path)
    bonds, dihedrals = topology(tris)
    a, v = area_volume(verts, tris)
    if len(tris) != 2 * len(verts) - 4 or len(bonds) != 3 * len(verts) - 6:
        raise SystemExit("mesh is not a closed surface of genus 0")
    if abs(a - area) > tol * area or abs(v - volume) > tol * volume:
        raise SystemExit("mesh area %.4f / volume %.4f differ from %.4f / %.4f"
                         % (a, v, area, volume))
    return a, v


def write_walls(path, blocks, velocities=None, group="mobile", fix_id="walls", slip=None,
                round=None):
    """Write the fix wall/bounceback command that matches the wall blocks.

    slip names a group (for example the membrane) that is reflected like
    a mirror and slides along the walls.  With round, that group sees the
    edges of the blocks rounded with this radius."""
    def fmt(v):
        if v >= INF:
            return "INF"
        if v <= -INF:
            return "-INF"
        return "%g" % v
    lines = ["fix %s %s wall/bounceback" % (fix_id, group)]
    if slip:
        lines[0] += " slip %s" % slip
    if round:
        lines[0] += " round %g" % round
    for i, b in enumerate(blocks):
        words = ["block"] + [fmt(v) for v in b]
        if velocities is not None and velocities[i] is not None:
            words += ["velocity"] + ["%g" % u for u in velocities[i]]
        lines.append("    " + " ".join(words))
    with open(path, "w") as f:
        f.write("# written by make_data.py: wall blocks of this box\n")
        f.write(" &\n".join(lines) + "\n")


