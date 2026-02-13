"""
MagSafe Mountain Stand v2 - STL Exporter (zero dependencies)
Fixed: wider/rounder mountain, flush MagSafe pocket, no floating geometry,
       integrated rock texture, proper proportions matching concept render.
Units: mm
"""
import struct
import math
import os
import random

# =============================================
# Parameters
# =============================================
BASE_D = 130.0
BASE_R = BASE_D / 2
BASE_H = 8.0
BASE_CHAMFER = 2.0

# Mountain: wider, rounder, squatter
MTN_H = 75.0           # shorter than v1
MTN_RX = 55.0           # wider
MTN_RY = 45.0           # deeper
MTN_PEAK_OFFSET_X = -3.0
MTN_PEAK_OFFSET_Y = 5.0  # peak toward the back
MTN_Y_CENTER = 8.0       # mountain center offset toward back

# MagSafe pocket - recessed into front face
MS_POCKET_D = 57.0
MS_POCKET_R = MS_POCKET_D / 2
MS_POCKET_DEPTH = 7.0
MS_CENTER_H = 45.0       # lower on the face
MS_TILT_DEG = 75.0       # angle from horizontal

# Cable
CABLE_D = 5.0

# Cross - directly on peak, no sphere
CROSS_H = 22.0
CROSS_ARM_W = 14.0
CROSS_BAR = 2.8
CROSS_ARM_POS = 0.68     # arm at 68% of cross height

# Coin tray - recessed into base surface
TRAY_W = 28.0
TRAY_D = 18.0
TRAY_DEPTH = 5.0
TRAY_X = -32.0
TRAY_Y = -18.0

# Book
BOOK_W = 32.0
BOOK_D = 22.0
BOOK_H = 3.0
BOOK_X = 12.0
BOOK_Y = -32.0

SEG = 48
random.seed(42)  # deterministic rock texture

# =============================================
# Geometry helpers
# =============================================

def compute_normal(tri):
    (x0,y0,z0),(x1,y1,z1),(x2,y2,z2) = tri
    ux,uy,uz = x1-x0,y1-y0,z1-z0
    vx,vy,vz = x2-x0,y2-y0,z2-z0
    nx = uy*vz-uz*vy; ny = uz*vx-ux*vz; nz = ux*vy-uy*vx
    l = math.sqrt(nx*nx+ny*ny+nz*nz)
    return (nx/l,ny/l,nz/l) if l>0 else (0,0,1)

def write_binary_stl(filename, triangles):
    with open(filename,'wb') as f:
        f.write(b'\x00'*80)
        f.write(struct.pack('<I',len(triangles)))
        for tri in triangles:
            n = compute_normal(tri)
            f.write(struct.pack('<fff',*n))
            for v in tri: f.write(struct.pack('<fff',*v))
            f.write(struct.pack('<H',0))

def box_tris(x,y,z,dx,dy,dz):
    v = [(x,y,z),(x+dx,y,z),(x+dx,y+dy,z),(x,y+dy,z),
         (x,y,z+dz),(x+dx,y,z+dz),(x+dx,y+dy,z+dz),(x,y+dy,z+dz)]
    faces = [(0,2,1),(0,3,2),(4,5,6),(4,6,7),(0,1,5),(0,5,4),
             (2,3,7),(2,7,6),(0,4,7),(0,7,3),(1,2,6),(1,6,5)]
    return [(v[a],v[b],v[c]) for a,b,c in faces]


# =============================================
# Mountain surface with integrated rock texture
# =============================================

def mountain_profile(t):
    """Mountain radius factor at height t (0=base, 1=peak).
    Uses a smooth curve that's wide at the bottom and rounds to a peak.
    """
    # Combination of cosine and power for natural mountain shape
    return max(0, math.cos(t * math.pi / 2) ** 0.6)


def rock_noise(angle, height_t, seed_offset=0):
    """Procedural rock noise integrated into the surface.
    Returns a multiplier (0.85 to 1.15) for radius perturbation.
    """
    # Multiple frequencies for realistic rock texture
    n = 0.0
    n += 0.08 * math.sin(angle * 7 + seed_offset + height_t * 12)
    n += 0.05 * math.sin(angle * 13 + seed_offset * 2.3 + height_t * 8)
    n += 0.04 * math.sin(angle * 19 + seed_offset * 0.7 + height_t * 20)
    # Larger rock ledge features
    n += 0.06 * max(0, math.sin(angle * 3 + height_t * 5)) ** 2
    # Reduce noise at very top and very bottom for cleaner transitions
    envelope = math.sin(max(0.05, min(0.95, height_t)) * math.pi)
    return 1.0 + n * envelope


def build_mountain_surface():
    """Build mountain as stacked elliptical slices with rock texture."""
    tris = []
    n_slices = 30

    for s in range(n_slices):
        t0 = s / n_slices
        t1 = (s + 1) / n_slices
        z0 = BASE_H + t0 * MTN_H
        z1 = BASE_H + t1 * MTN_H

        prof0 = mountain_profile(t0)
        prof1 = mountain_profile(t1)

        # Peak offset interpolation
        ox0 = MTN_PEAK_OFFSET_X * t0
        oy0 = MTN_PEAK_OFFSET_Y * t0
        ox1 = MTN_PEAK_OFFSET_X * t1
        oy1 = MTN_PEAK_OFFSET_Y * t1

        for i in range(SEG):
            a0 = 2 * math.pi * i / SEG
            a1 = 2 * math.pi * (i + 1) / SEG

            # Rock noise at each vertex
            rn00 = rock_noise(a0, t0, 0)
            rn01 = rock_noise(a1, t0, 0)
            rn10 = rock_noise(a0, t1, 0)
            rn11 = rock_noise(a1, t1, 0)

            # Vertex positions
            rx0 = MTN_RX * prof0
            ry0 = MTN_RY * prof0
            rx1 = MTN_RX * prof1
            ry1 = MTN_RY * prof1

            b0 = (ox0 + rx0*rn00*math.cos(a0), MTN_Y_CENTER + oy0 + ry0*rn00*math.sin(a0), z0)
            b1 = (ox0 + rx0*rn01*math.cos(a1), MTN_Y_CENTER + oy0 + ry0*rn01*math.sin(a1), z0)
            tp0 = (ox1 + rx1*rn10*math.cos(a0), MTN_Y_CENTER + oy1 + ry1*rn10*math.sin(a0), z1)
            tp1 = (ox1 + rx1*rn11*math.cos(a1), MTN_Y_CENTER + oy1 + ry1*rn11*math.sin(a1), z1)

            tris.append((b0, b1, tp1))
            tris.append((b0, tp1, tp0))

    # Top cap (small disk to close the peak)
    peak_z = BASE_H + MTN_H
    peak_r = 4.0  # small rounded peak
    peak_cx = MTN_PEAK_OFFSET_X
    peak_cy = MTN_Y_CENTER + MTN_PEAK_OFFSET_Y
    for i in range(SEG):
        a0 = 2*math.pi*i/SEG
        a1 = 2*math.pi*(i+1)/SEG
        c = (peak_cx, peak_cy, peak_z)
        p0 = (peak_cx + peak_r*math.cos(a0), peak_cy + peak_r*math.sin(a0), peak_z)
        p1 = (peak_cx + peak_r*math.cos(a1), peak_cy + peak_r*math.sin(a1), peak_z)
        tris.append((c, p0, p1))

    return tris


# =============================================
# MagSafe pocket (boolean-style recess)
# =============================================

def build_magsafe_pocket():
    """Build the MagSafe pocket as geometry that visually represents the recess.
    Creates: pocket walls (cylinder into the mountain) + bottom disk.
    """
    tris = []
    tilt = math.radians(MS_TILT_DEG)

    # Pocket center on front face
    cx = 0
    cy = MTN_Y_CENTER - MTN_RY * 0.55  # on the front face
    cz = BASE_H + MS_CENTER_H

    # Direction into mountain (tilted back)
    dy_dir = math.cos(tilt)
    dz_dir = -math.sin(tilt)

    # Build pocket as a cylinder: entrance ring -> bottom
    for i in range(SEG):
        a0 = 2*math.pi*i/SEG
        a1 = 2*math.pi*(i+1)/SEG

        r = MS_POCKET_R

        # Entrance points (on mountain surface)
        # The pocket circle is perpendicular to the tilt direction
        # "up" in pocket space = along tilt normal
        e0 = (
            cx + r * math.cos(a0),
            cy + r * math.sin(a0) * math.sin(tilt),
            cz + r * math.sin(a0) * math.cos(tilt)
        )
        e1 = (
            cx + r * math.cos(a1),
            cy + r * math.sin(a1) * math.sin(tilt),
            cz + r * math.sin(a1) * math.cos(tilt)
        )

        # Bottom points (pushed into mountain)
        d = MS_POCKET_DEPTH
        b0 = (e0[0], e0[1] + d*dy_dir, e0[2] + d*dz_dir)
        b1 = (e1[0], e1[1] + d*dy_dir, e1[2] + d*dz_dir)

        # Pocket inner walls (normals pointing inward = visible from inside)
        tris.append((e0, b0, b1))
        tris.append((e0, b1, e1))

        # Bottom disk
        bottom_center = (cx, cy + d*dy_dir, cz + d*dz_dir)
        tris.append((bottom_center, b1, b0))

    # Entrance rim ring (to seal the pocket opening cleanly)
    # Slightly larger ring to create a gold rim like in concept
    rim_width = 3.0
    for i in range(SEG):
        a0 = 2*math.pi*i/SEG
        a1 = 2*math.pi*(i+1)/SEG

        ri = MS_POCKET_R
        ro = MS_POCKET_R + rim_width

        # Inner ring points
        i0 = (cx + ri*math.cos(a0), cy + ri*math.sin(a0)*math.sin(tilt), cz + ri*math.sin(a0)*math.cos(tilt))
        i1 = (cx + ri*math.cos(a1), cy + ri*math.sin(a1)*math.sin(tilt), cz + ri*math.sin(a1)*math.cos(tilt))
        # Outer ring points
        o0 = (cx + ro*math.cos(a0), cy + ro*math.sin(a0)*math.sin(tilt), cz + ro*math.sin(a0)*math.cos(tilt))
        o1 = (cx + ro*math.cos(a1), cy + ro*math.sin(a1)*math.sin(tilt), cz + ro*math.sin(a1)*math.cos(tilt))

        # Ring face (pointing outward)
        tris.append((i0, o0, o1))
        tris.append((i0, o1, i1))

    return tris


# =============================================
# Coin tray (recessed into base)
# =============================================

def build_coin_tray():
    """Coin tray as a recessed pocket in the base surface."""
    tris = []
    x, y, z = TRAY_X, TRAY_Y, BASE_H
    w, d, h = TRAY_W, TRAY_D, TRAY_DEPTH

    # Floor of tray
    tris.append(((x,y,z-h+2), (x+w,y,z-h+2), (x+w,y+d,z-h+2)))
    tris.append(((x,y,z-h+2), (x+w,y+d,z-h+2), (x,y+d,z-h+2)))

    # Inner walls (visible from above)
    # Front wall
    tris.append(((x,y,z-h+2), (x,y,z), (x+w,y,z)))
    tris.append(((x,y,z-h+2), (x+w,y,z), (x+w,y,z-h+2)))
    # Back wall
    tris.append(((x,y+d,z-h+2), (x+w,y+d,z), (x,y+d,z)))
    tris.append(((x,y+d,z-h+2), (x+w,y+d,z-h+2), (x+w,y+d,z)))
    # Left wall
    tris.append(((x,y,z-h+2), (x,y+d,z), (x,y,z)))
    tris.append(((x,y,z-h+2), (x,y+d,z-h+2), (x,y+d,z)))
    # Right wall
    tris.append(((x+w,y,z-h+2), (x+w,y,z), (x+w,y+d,z)))
    tris.append(((x+w,y,z-h+2), (x+w,y+d,z), (x+w,y+d,z-h+2)))

    # Thin raised rim around the tray
    rim = 1.2
    rim_h = 1.5
    tris += box_tris(x-rim, y-rim, z, w+2*rim, rim, rim_h)
    tris += box_tris(x-rim, y+d, z, w+2*rim, rim, rim_h)
    tris += box_tris(x-rim, y, z, rim, d, rim_h)
    tris += box_tris(x+w, y, z, rim, d, rim_h)

    return tris


# =============================================
# Open book
# =============================================

def build_book():
    """Open book sitting on the base."""
    tris = []
    bx, by, bz = BOOK_X, BOOK_Y, BASE_H

    # Left page (slight upward angle)
    hw = BOOK_W / 2 - 1
    # Left page - flat for simplicity
    tris += box_tris(bx - BOOK_W/2, by, bz, hw, BOOK_D, BOOK_H)
    # Right page
    tris += box_tris(bx + 1, by, bz, hw, BOOK_D, BOOK_H)
    # Spine (slightly raised ridge)
    tris += box_tris(bx - 1, by, bz, 2, BOOK_D, BOOK_H + 0.8)

    return tris


# =============================================
# Cross (directly on peak)
# =============================================

def build_cross():
    """Cross sitting directly on the mountain peak."""
    tris = []
    cx = MTN_PEAK_OFFSET_X
    cy = MTN_Y_CENTER + MTN_PEAK_OFFSET_Y
    cz = BASE_H + MTN_H  # directly on peak

    arm_h = CROSS_H * CROSS_ARM_POS

    # Vertical bar
    tris += box_tris(cx-CROSS_BAR/2, cy-CROSS_BAR/2, cz, CROSS_BAR, CROSS_BAR, CROSS_H)
    # Horizontal arm
    tris += box_tris(cx-CROSS_ARM_W/2, cy-CROSS_BAR/2, cz+arm_h, CROSS_ARM_W, CROSS_BAR, CROSS_BAR)

    return tris


# =============================================
# Base platform
# =============================================

def build_base():
    """Chamfered circular base."""
    tris = []
    # Bottom face (smaller for chamfer)
    r_bottom = BASE_R - BASE_CHAMFER
    for i in range(SEG):
        a0 = 2*math.pi*i/SEG
        a1 = 2*math.pi*(i+1)/SEG
        c = (0, 0, 0)
        p0 = (r_bottom*math.cos(a0), r_bottom*math.sin(a0), 0)
        p1 = (r_bottom*math.cos(a1), r_bottom*math.sin(a1), 0)
        tris.append((c, p1, p0))

    # Chamfer (angled ring from bottom to main body)
    for i in range(SEG):
        a0 = 2*math.pi*i/SEG
        a1 = 2*math.pi*(i+1)/SEG
        b0 = (r_bottom*math.cos(a0), r_bottom*math.sin(a0), 0)
        b1 = (r_bottom*math.cos(a1), r_bottom*math.sin(a1), 0)
        t0 = (BASE_R*math.cos(a0), BASE_R*math.sin(a0), BASE_CHAMFER)
        t1 = (BASE_R*math.cos(a1), BASE_R*math.sin(a1), BASE_CHAMFER)
        tris.append((b0, b1, t1))
        tris.append((b0, t1, t0))

    # Side wall
    for i in range(SEG):
        a0 = 2*math.pi*i/SEG
        a1 = 2*math.pi*(i+1)/SEG
        b0 = (BASE_R*math.cos(a0), BASE_R*math.sin(a0), BASE_CHAMFER)
        b1 = (BASE_R*math.cos(a1), BASE_R*math.sin(a1), BASE_CHAMFER)
        t0 = (BASE_R*math.cos(a0), BASE_R*math.sin(a0), BASE_H)
        t1 = (BASE_R*math.cos(a1), BASE_R*math.sin(a1), BASE_H)
        tris.append((b0, b1, t1))
        tris.append((b0, t1, t0))

    # Top face
    for i in range(SEG):
        a0 = 2*math.pi*i/SEG
        a1 = 2*math.pi*(i+1)/SEG
        c = (0, 0, BASE_H)
        p0 = (BASE_R*math.cos(a0), BASE_R*math.sin(a0), BASE_H)
        p1 = (BASE_R*math.cos(a1), BASE_R*math.sin(a1), BASE_H)
        tris.append((c, p0, p1))

    return tris


# =============================================
# Assemble everything
# =============================================

def build_stand():
    tris = []
    tris += build_base()
    tris += build_mountain_surface()
    tris += build_magsafe_pocket()
    tris += build_cross()
    tris += build_coin_tray()
    tris += build_book()
    return tris


if __name__ == "__main__":
    print("Building MagSafe Mountain Stand v2...")
    total_h = BASE_H + MTN_H + CROSS_H
    print(f"  Size: {BASE_D:.0f}mm dia x {total_h:.0f}mm tall")
    print(f"  AD5X check: {'PASS' if max(BASE_D, total_h) <= 220 else 'FAIL'}")

    tris = build_stand()
    print(f"  {len(tris)} triangles")

    out_dir = os.path.dirname(os.path.abspath(__file__))
    stl_path = os.path.join(out_dir, "magsafe_mountain_stand.stl")
    write_binary_stl(stl_path, tris)

    size_kb = os.path.getsize(stl_path) / 1024
    print(f"  Exported: {stl_path}")
    print(f"  Size: {size_kb:.1f} KB")
    print()
    print("Changes from v1:")
    print("  - Mountain: wider, rounder, squatter profile")
    print("  - MagSafe: flush pocket recess with rim ring")
    print("  - Rock texture: integrated into surface (no floating spheres)")
    print("  - Cross: sits directly on peak (no sphere)")
    print("  - Coin tray: recessed into base with thin rim")
    print("  - Book: flush on base surface")
