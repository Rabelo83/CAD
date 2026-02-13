"""
LEGO 2x4 Brick - STL Exporter (zero dependencies)
Generates a binary STL file from parametric dimensions.
Units: mm
"""
import struct
import math
import os

# === Parameters (match the .scad file) ===
STUD_PITCH = 8.0
STUDS_X = 4
STUDS_Y = 2

BRICK_L = STUDS_X * STUD_PITCH - 0.2  # 31.8
BRICK_W = STUDS_Y * STUD_PITCH - 0.2  # 15.8
BRICK_H = 11.4

WALL = 1.5
FLOOR_T = 1.0

STUD_D = 4.8
STUD_H = 1.8
STUD_R = STUD_D / 2

TUBE_OD = 6.51
TUBE_ID = 4.8
TUBE_OR = TUBE_OD / 2
TUBE_IR = TUBE_ID / 2

CYL_SEGMENTS = 48  # smoothness for cylinders


def cylinder_triangles(cx, cy, z_bot, z_top, radius, segments, outward=True):
    """Generate triangles for a cylinder wall."""
    tris = []
    for i in range(segments):
        a0 = 2 * math.pi * i / segments
        a1 = 2 * math.pi * (i + 1) / segments
        x0, y0 = cx + radius * math.cos(a0), cy + radius * math.sin(a0)
        x1, y1 = cx + radius * math.cos(a1), cy + radius * math.sin(a1)

        if outward:
            tris.append(((x0, y0, z_bot), (x1, y1, z_bot), (x1, y1, z_top)))
            tris.append(((x0, y0, z_bot), (x1, y1, z_top), (x0, y0, z_top)))
        else:
            tris.append(((x0, y0, z_bot), (x1, y1, z_top), (x1, y1, z_bot)))
            tris.append(((x0, y0, z_bot), (x0, y0, z_top), (x1, y1, z_top)))
    return tris


def disk_triangles(cx, cy, z, radius, segments, normal_up=True):
    """Generate triangles for a circular disk (top/bottom of cylinder)."""
    tris = []
    for i in range(segments):
        a0 = 2 * math.pi * i / segments
        a1 = 2 * math.pi * (i + 1) / segments
        x0, y0 = cx + radius * math.cos(a0), cy + radius * math.sin(a0)
        x1, y1 = cx + radius * math.cos(a1), cy + radius * math.sin(a1)
        if normal_up:
            tris.append(((cx, cy, z), (x0, y0, z), (x1, y1, z)))
        else:
            tris.append(((cx, cy, z), (x1, y1, z), (x0, y0, z)))
    return tris


def annular_ring_triangles(cx, cy, z, r_inner, r_outer, segments, normal_up=True):
    """Generate triangles for a ring (annulus) face."""
    tris = []
    for i in range(segments):
        a0 = 2 * math.pi * i / segments
        a1 = 2 * math.pi * (i + 1) / segments
        ox0 = cx + r_outer * math.cos(a0)
        oy0 = cy + r_outer * math.sin(a0)
        ox1 = cx + r_outer * math.cos(a1)
        oy1 = cy + r_outer * math.sin(a1)
        ix0 = cx + r_inner * math.cos(a0)
        iy0 = cy + r_inner * math.sin(a0)
        ix1 = cx + r_inner * math.cos(a1)
        iy1 = cy + r_inner * math.sin(a1)
        if normal_up:
            tris.append(((ix0, iy0, z), (ox0, oy0, z), (ox1, oy1, z)))
            tris.append(((ix0, iy0, z), (ox1, oy1, z), (ix1, iy1, z)))
        else:
            tris.append(((ix0, iy0, z), (ox1, oy1, z), (ox0, oy0, z)))
            tris.append(((ix0, iy0, z), (ix1, iy1, z), (ox1, oy1, z)))
    return tris


def box_triangles(x, y, z, dx, dy, dz):
    """Generate 12 triangles for an axis-aligned box."""
    v = [
        (x, y, z), (x+dx, y, z), (x+dx, y+dy, z), (x, y+dy, z),           # bottom
        (x, y, z+dz), (x+dx, y, z+dz), (x+dx, y+dy, z+dz), (x, y+dy, z+dz)  # top
    ]
    faces = [
        (0,3,2), (0,2,1),  # bottom (z-)
        (4,5,6), (4,6,7),  # top (z+)
        (0,1,5), (0,5,4),  # front (y-)
        (2,3,7), (2,7,6),  # back (y+)
        (0,4,7), (0,7,3),  # left (x-)
        (1,2,6), (1,6,5),  # right (x+)
    ]
    return [(v[a], v[b], v[c]) for a, b, c in faces]


def compute_normal(tri):
    """Compute face normal from triangle vertices."""
    (x0,y0,z0), (x1,y1,z1), (x2,y2,z2) = tri
    ux, uy, uz = x1-x0, y1-y0, z1-z0
    vx, vy, vz = x2-x0, y2-y0, z2-z0
    nx = uy*vz - uz*vy
    ny = uz*vx - ux*vz
    nz = ux*vy - uy*vx
    length = math.sqrt(nx*nx + ny*ny + nz*nz)
    if length > 0:
        nx, ny, nz = nx/length, ny/length, nz/length
    return (nx, ny, nz)


def write_binary_stl(filename, triangles):
    """Write triangles to a binary STL file."""
    with open(filename, 'wb') as f:
        f.write(b'\x00' * 80)  # header
        f.write(struct.pack('<I', len(triangles)))
        for tri in triangles:
            nx, ny, nz = compute_normal(tri)
            f.write(struct.pack('<fff', nx, ny, nz))
            for vertex in tri:
                f.write(struct.pack('<fff', *vertex))
            f.write(struct.pack('<H', 0))  # attribute byte count


def build_lego_brick():
    """Build all triangles for the LEGO 2x4 brick."""
    triangles = []

    # --- Outer shell (hollow box) ---
    # Outer box
    triangles += box_triangles(0, 0, 0, BRICK_L, BRICK_W, BRICK_H)
    # Inner cavity (inverted normals = subtract)
    inner = box_triangles(WALL, WALL, FLOOR_T,
                          BRICK_L - 2*WALL, BRICK_W - 2*WALL, BRICK_H - FLOOR_T + 0.01)
    # Flip winding to make inner cavity
    triangles += [(c, b, a) for a, b, c in inner]

    # Top face: cut out the inner rectangle by replacing the solid top
    # The box_triangles already created a solid top, but we need to hollow it.
    # For simplicity with pure triangles, we'll create the shell as walls + floor + top rim.

    # Let's rebuild more carefully:
    triangles = []

    # Bottom face (solid)
    triangles += [
        ((0,0,0), (BRICK_L,BRICK_W,0), (BRICK_L,0,0)),
        ((0,0,0), (0,BRICK_W,0), (BRICK_L,BRICK_W,0)),
    ]

    # Top face (rim only - outer minus inner)
    # Outer top rectangle edges -> inner top rectangle edges
    ox0, oy0 = 0, 0
    ox1, oy1 = BRICK_L, BRICK_W
    ix0, iy0 = WALL, WALL
    ix1, iy1 = BRICK_L - WALL, BRICK_W - WALL
    zt = BRICK_H

    # Front rim strip
    triangles += [
        ((ox0,oy0,zt), (ox1,oy0,zt), (ix1,iy0,zt)),
        ((ox0,oy0,zt), (ix1,iy0,zt), (ix0,iy0,zt)),
    ]
    # Back rim strip
    triangles += [
        ((ox0,oy1,zt), (ix0,iy1,zt), (ix1,iy1,zt)),
        ((ox0,oy1,zt), (ix1,iy1,zt), (ox1,oy1,zt)),
    ]
    # Left rim strip
    triangles += [
        ((ox0,oy0,zt), (ix0,iy0,zt), (ix0,iy1,zt)),
        ((ox0,oy0,zt), (ix0,iy1,zt), (ox0,oy1,zt)),
    ]
    # Right rim strip
    triangles += [
        ((ox1,oy0,zt), (ox1,oy1,zt), (ix1,iy1,zt)),
        ((ox1,oy0,zt), (ix1,iy1,zt), (ix1,iy0,zt)),
    ]

    # Outer walls (4 sides)
    # Front wall (y=0)
    triangles += [
        ((0,0,0), (BRICK_L,0,0), (BRICK_L,0,BRICK_H)),
        ((0,0,0), (BRICK_L,0,BRICK_H), (0,0,BRICK_H)),
    ]
    # Back wall (y=W)
    triangles += [
        ((0,BRICK_W,0), (0,BRICK_W,BRICK_H), (BRICK_L,BRICK_W,BRICK_H)),
        ((0,BRICK_W,0), (BRICK_L,BRICK_W,BRICK_H), (BRICK_L,BRICK_W,0)),
    ]
    # Left wall (x=0)
    triangles += [
        ((0,0,0), (0,0,BRICK_H), (0,BRICK_W,BRICK_H)),
        ((0,0,0), (0,BRICK_W,BRICK_H), (0,BRICK_W,0)),
    ]
    # Right wall (x=L)
    triangles += [
        ((BRICK_L,0,0), (BRICK_L,BRICK_W,0), (BRICK_L,BRICK_W,BRICK_H)),
        ((BRICK_L,0,0), (BRICK_L,BRICK_W,BRICK_H), (BRICK_L,0,BRICK_H)),
    ]

    # Inner walls (4 sides of cavity)
    zb = FLOOR_T
    # Front inner (y=WALL)
    triangles += [
        ((WALL,WALL,zb), (BRICK_L-WALL,WALL,BRICK_H), (BRICK_L-WALL,WALL,zb)),
        ((WALL,WALL,zb), (WALL,WALL,BRICK_H), (BRICK_L-WALL,WALL,BRICK_H)),
    ]
    # Back inner (y=W-WALL)
    bw = BRICK_W - WALL
    triangles += [
        ((WALL,bw,zb), (BRICK_L-WALL,bw,zb), (BRICK_L-WALL,bw,BRICK_H)),
        ((WALL,bw,zb), (BRICK_L-WALL,bw,BRICK_H), (WALL,bw,BRICK_H)),
    ]
    # Left inner (x=WALL)
    triangles += [
        ((WALL,WALL,zb), (WALL,bw,BRICK_H), (WALL,WALL,BRICK_H)),
        ((WALL,WALL,zb), (WALL,bw,zb), (WALL,bw,BRICK_H)),
    ]
    # Right inner (x=L-WALL)
    rw = BRICK_L - WALL
    triangles += [
        ((rw,WALL,zb), (rw,WALL,BRICK_H), (rw,bw,BRICK_H)),
        ((rw,WALL,zb), (rw,bw,BRICK_H), (rw,bw,zb)),
    ]

    # Inner floor (ceiling of cavity)
    triangles += [
        ((WALL,WALL,zb), (BRICK_L-WALL,WALL,zb), (BRICK_L-WALL,bw,zb)),
        ((WALL,WALL,zb), (BRICK_L-WALL,bw,zb), (WALL,bw,zb)),
    ]

    # --- Studs (solid cylinders on top) ---
    for ix in range(STUDS_X):
        for iy in range(STUDS_Y):
            cx = STUD_PITCH/2 + ix * STUD_PITCH
            cy = STUD_PITCH/2 + iy * STUD_PITCH
            z_bot = BRICK_H
            z_top = BRICK_H + STUD_H

            # Cylinder wall
            triangles += cylinder_triangles(cx, cy, z_bot, z_top, STUD_R, CYL_SEGMENTS, outward=True)
            # Top cap
            triangles += disk_triangles(cx, cy, z_top, STUD_R, CYL_SEGMENTS, normal_up=True)
            # Bottom cap (where stud meets brick top - fills the hole)
            triangles += disk_triangles(cx, cy, z_bot, STUD_R, CYL_SEGMENTS, normal_up=False)

    # --- Underside tubes (hollow cylinders) ---
    # Tubes go between studs: (STUDS_X-1) x (STUDS_Y-1) grid
    for ix in range(STUDS_X - 1):
        for iy in range(STUDS_Y - 1):
            cx = STUD_PITCH + ix * STUD_PITCH
            cy = STUD_PITCH + iy * STUD_PITCH
            z_bot = FLOOR_T
            z_top = BRICK_H

            # Outer wall
            triangles += cylinder_triangles(cx, cy, z_bot, z_top, TUBE_OR, CYL_SEGMENTS, outward=True)
            # Inner wall (inward normals)
            triangles += cylinder_triangles(cx, cy, z_bot, z_top, TUBE_IR, CYL_SEGMENTS, outward=False)
            # Top ring
            triangles += annular_ring_triangles(cx, cy, z_top, TUBE_IR, TUBE_OR, CYL_SEGMENTS, normal_up=True)
            # Bottom ring
            triangles += annular_ring_triangles(cx, cy, z_bot, TUBE_IR, TUBE_OR, CYL_SEGMENTS, normal_up=False)

    return triangles


if __name__ == "__main__":
    print("Building LEGO 2x4 brick geometry...")
    tris = build_lego_brick()
    print(f"  {len(tris)} triangles generated")

    out_dir = os.path.dirname(os.path.abspath(__file__))
    stl_path = os.path.join(out_dir, "lego_2x4_brick.stl")

    write_binary_stl(stl_path, tris)
    size_kb = os.path.getsize(stl_path) / 1024
    print(f"  Exported: {stl_path}")
    print(f"  File size: {size_kb:.1f} KB")
    print(f"\nOpen this STL in any online viewer:")
    print(f"  - Onshape (import)")
    print(f"  - FreeCAD")
    print(f"  - Any slicer (Cura, PrusaSlicer, etc.)")
