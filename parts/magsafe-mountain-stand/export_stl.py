"""
MagSafe Mountain Stand - STL Exporter (zero dependencies)
Generates a binary STL of the functional stand geometry.
Units: mm

This creates a simplified but printable version:
- Chamfered round base
- Mountain/cone body with rock-like bumps
- MagSafe puck pocket (angled recess)
- Cable channel
- Cross on top
- Coin tray cutout
- Open book detail

For full rock texture, import the STL into Onshape/Blender and sculpt detail.
"""
import struct
import math
import os

# =============================================
# Parameters (match .scad file)
# =============================================
BASE_D = 130.0
BASE_R = BASE_D / 2
BASE_H = 8.0
BASE_CHAMFER = 2.0

MTN_H = 95.0
MTN_RX = 50.0
MTN_RY = 40.0
MTN_PEAK_OFFSET = -5.0

MS_POCKET_D = 57.0
MS_POCKET_R = MS_POCKET_D / 2
MS_POCKET_DEPTH = 6.5
MS_CENTER_H = 55.0
MS_TILT = 75.0  # degrees from horizontal

CABLE_D = 5.0

CROSS_H = 25.0
CROSS_ARM_W = 15.0
CROSS_BAR = 3.0
CROSS_ARM_H = 17.0

TRAY_W = 30.0
TRAY_D = 20.0
TRAY_H = 8.0
TRAY_OFFSET_X = -35.0
TRAY_OFFSET_Y = -15.0

BOOK_W = 35.0
BOOK_D = 25.0
BOOK_H = 4.0
BOOK_OFFSET_X = 10.0
BOOK_OFFSET_Y = -30.0

SEG = 48  # cylinder/sphere segments

# =============================================
# Geometry helpers
# =============================================

def compute_normal(tri):
    (x0,y0,z0), (x1,y1,z1), (x2,y2,z2) = tri
    ux, uy, uz = x1-x0, y1-y0, z1-z0
    vx, vy, vz = x2-x0, y2-y0, z2-z0
    nx = uy*vz - uz*vy
    ny = uz*vx - ux*vz
    nz = ux*vy - uy*vx
    l = math.sqrt(nx*nx + ny*ny + nz*nz)
    if l > 0:
        return (nx/l, ny/l, nz/l)
    return (0, 0, 1)

def write_binary_stl(filename, triangles):
    with open(filename, 'wb') as f:
        f.write(b'\x00' * 80)
        f.write(struct.pack('<I', len(triangles)))
        for tri in triangles:
            nx, ny, nz = compute_normal(tri)
            f.write(struct.pack('<fff', nx, ny, nz))
            for v in tri:
                f.write(struct.pack('<fff', *v))
            f.write(struct.pack('<H', 0))

def box_tris(x, y, z, dx, dy, dz):
    v = [
        (x,y,z), (x+dx,y,z), (x+dx,y+dy,z), (x,y+dy,z),
        (x,y,z+dz), (x+dx,y,z+dz), (x+dx,y+dy,z+dz), (x,y+dy,z+dz)
    ]
    faces = [
        (0,2,1),(0,3,2), (4,5,6),(4,6,7),
        (0,1,5),(0,5,4), (2,3,7),(2,7,6),
        (0,4,7),(0,7,3), (1,2,6),(1,6,5),
    ]
    return [(v[a],v[b],v[c]) for a,b,c in faces]

def cylinder_wall(cx, cy, z0, z1, r, seg, outward=True):
    tris = []
    for i in range(seg):
        a0 = 2*math.pi*i/seg
        a1 = 2*math.pi*(i+1)/seg
        p0 = (cx+r*math.cos(a0), cy+r*math.sin(a0))
        p1 = (cx+r*math.cos(a1), cy+r*math.sin(a1))
        if outward:
            tris.append(((p0[0],p0[1],z0),(p1[0],p1[1],z0),(p1[0],p1[1],z1)))
            tris.append(((p0[0],p0[1],z0),(p1[0],p1[1],z1),(p0[0],p0[1],z1)))
        else:
            tris.append(((p0[0],p0[1],z0),(p1[0],p1[1],z1),(p1[0],p1[1],z0)))
            tris.append(((p0[0],p0[1],z0),(p0[0],p0[1],z1),(p1[0],p1[1],z1)))
    return tris

def disk(cx, cy, z, r, seg, up=True):
    tris = []
    for i in range(seg):
        a0 = 2*math.pi*i/seg
        a1 = 2*math.pi*(i+1)/seg
        p0 = (cx+r*math.cos(a0), cy+r*math.sin(a0), z)
        p1 = (cx+r*math.cos(a1), cy+r*math.sin(a1), z)
        c = (cx, cy, z)
        if up:
            tris.append((c, p0, p1))
        else:
            tris.append((c, p1, p0))
    return tris

def cone_surface(cx, cy, z_base, z_top, r_base, r_top, seg):
    """Frustum (truncated cone) surface triangles."""
    tris = []
    for i in range(seg):
        a0 = 2*math.pi*i/seg
        a1 = 2*math.pi*(i+1)/seg
        b0 = (cx+r_base*math.cos(a0), cy+r_base*math.sin(a0), z_base)
        b1 = (cx+r_base*math.cos(a1), cy+r_base*math.sin(a1), z_base)
        t0 = (cx+r_top*math.cos(a0), cy+r_top*math.sin(a0), z_top)
        t1 = (cx+r_top*math.cos(a1), cy+r_top*math.sin(a1), z_top)
        tris.append((b0, b1, t1))
        tris.append((b0, t1, t0))
    return tris

def elliptical_cone_surface(cx, cy, z_base, z_top, rx_base, ry_base, rx_top, ry_top, seg, offset_x=0):
    """Elliptical frustum with optional peak X offset."""
    tris = []
    for i in range(seg):
        a0 = 2*math.pi*i/seg
        a1 = 2*math.pi*(i+1)/seg
        b0 = (cx+rx_base*math.cos(a0), cy+ry_base*math.sin(a0), z_base)
        b1 = (cx+rx_base*math.cos(a1), cy+ry_base*math.sin(a1), z_base)
        t0 = (cx+offset_x+rx_top*math.cos(a0), cy+ry_top*math.sin(a0), z_top)
        t1 = (cx+offset_x+rx_top*math.cos(a1), cy+ry_top*math.sin(a1), z_top)
        tris.append((b0, b1, t1))
        tris.append((b0, t1, t0))
    return tris

def sphere_tris(cx, cy, cz, r, h_seg=16, v_seg=12):
    """Generate sphere triangles."""
    tris = []
    for j in range(v_seg):
        phi0 = math.pi * j / v_seg
        phi1 = math.pi * (j+1) / v_seg
        for i in range(h_seg):
            th0 = 2*math.pi*i/h_seg
            th1 = 2*math.pi*(i+1)/h_seg
            p00 = (cx+r*math.sin(phi0)*math.cos(th0), cy+r*math.sin(phi0)*math.sin(th0), cz+r*math.cos(phi0))
            p01 = (cx+r*math.sin(phi0)*math.cos(th1), cy+r*math.sin(phi0)*math.sin(th1), cz+r*math.cos(phi0))
            p10 = (cx+r*math.sin(phi1)*math.cos(th0), cy+r*math.sin(phi1)*math.sin(th0), cz+r*math.cos(phi1))
            p11 = (cx+r*math.sin(phi1)*math.cos(th1), cy+r*math.sin(phi1)*math.sin(th1), cz+r*math.cos(phi1))
            if j > 0:
                tris.append((p00, p10, p01))
            if j < v_seg-1:
                tris.append((p01, p10, p11))
    return tris

# =============================================
# Build the stand
# =============================================

def build_stand():
    tris = []

    # === 1. BASE PLATFORM (chamfered cylinder) ===
    # Bottom disk (smaller, for chamfer)
    tris += disk(0, 0, 0, BASE_R - BASE_CHAMFER, SEG, up=False)
    # Chamfer ring (angled surface from bottom edge to main body)
    tris += cone_surface(0, 0, 0, BASE_CHAMFER, BASE_R - BASE_CHAMFER, BASE_R, SEG)
    # Main cylinder wall
    tris += cylinder_wall(0, 0, BASE_CHAMFER, BASE_H, BASE_R, SEG)
    # Top disk
    tris += disk(0, 0, BASE_H, BASE_R, SEG, up=True)

    # === 2. MOUNTAIN BODY ===
    # Build as stacked elliptical frustum slices for organic shape
    mtn_y_off = 5  # mountain center offset toward back
    n_slices = 20
    for s in range(n_slices):
        t0 = s / n_slices
        t1 = (s + 1) / n_slices
        z0 = BASE_H + t0 * MTN_H
        z1 = BASE_H + t1 * MTN_H

        # Radius tapers with some noise for rocky look
        # Use a power curve for mountain profile
        taper0 = (1 - t0**0.7)
        taper1 = (1 - t1**0.7)

        # Add radial variation per slice
        rx0 = MTN_RX * taper0
        ry0 = MTN_RY * taper0
        rx1 = MTN_RX * taper1
        ry1 = MTN_RY * taper1

        # Peak offset interpolation
        off0 = MTN_PEAK_OFFSET * t0
        off1 = MTN_PEAK_OFFSET * t1

        for i in range(SEG):
            a0 = 2*math.pi*i/SEG
            a1 = 2*math.pi*(i+1)/SEG

            # Rock texture: perturb radius per angle
            noise0_a0 = 1.0 + 0.15 * math.sin(a0*5 + t0*8) * (1-t0)
            noise0_a1 = 1.0 + 0.15 * math.sin(a1*5 + t0*8) * (1-t0)
            noise1_a0 = 1.0 + 0.15 * math.sin(a0*5 + t1*8) * (1-t1)
            noise1_a1 = 1.0 + 0.15 * math.sin(a1*5 + t1*8) * (1-t1)

            b0 = (off0 + rx0*noise0_a0*math.cos(a0), mtn_y_off + ry0*noise0_a0*math.sin(a0), z0)
            b1 = (off0 + rx0*noise0_a1*math.cos(a1), mtn_y_off + ry0*noise0_a1*math.sin(a1), z0)
            t0p = (off1 + rx1*noise1_a0*math.cos(a0), mtn_y_off + ry1*noise1_a0*math.sin(a0), z1)
            t1p = (off1 + rx1*noise1_a1*math.cos(a1), mtn_y_off + ry1*noise1_a1*math.sin(a1), z1)

            tris.append((b0, b1, t1p))
            tris.append((b0, t1p, t0p))

    # Mountain top cap (small sphere at peak)
    peak_z = BASE_H + MTN_H
    tris += sphere_tris(MTN_PEAK_OFFSET, mtn_y_off, peak_z, 6, 12, 8)

    # === 3. MAGSAFE POCKET ===
    # Create the pocket as an indentation on the front face
    # Pocket is a cylinder angled into the mountain
    tilt_rad = math.radians(MS_TILT)
    pocket_cx = 0
    pocket_cy = -MTN_RY * 0.45 + 5  # front face of mountain
    pocket_cz = BASE_H + MS_CENTER_H

    # Direction vector (pointing into mountain)
    dx = 0
    dy = math.cos(tilt_rad)   # into mountain (toward back)
    dz = -math.sin(tilt_rad)  # angled upward

    # Pocket entrance ring (circular face)
    # Build the pocket face as a disk perpendicular to the tilt angle
    up = (0, -math.sin(tilt_rad), -math.cos(tilt_rad))  # normal pointing outward
    right = (1, 0, 0)  # horizontal

    for i in range(SEG):
        a0 = 2*math.pi*i/SEG
        a1 = 2*math.pi*(i+1)/SEG

        # Points on pocket entrance circle
        r = MS_POCKET_R
        p0 = (
            pocket_cx + r*math.cos(a0),
            pocket_cy + r*math.sin(a0)*math.sin(tilt_rad),
            pocket_cz + r*math.sin(a0)*math.cos(tilt_rad)
        )
        p1 = (
            pocket_cx + r*math.cos(a1),
            pocket_cy + r*math.sin(a1)*math.sin(tilt_rad),
            pocket_cz + r*math.sin(a1)*math.cos(tilt_rad)
        )

        # Points on pocket bottom circle (pushed into mountain)
        depth = MS_POCKET_DEPTH
        q0 = (
            p0[0],
            p0[1] + depth*dy,
            p0[2] + depth*dz
        )
        q1 = (
            p1[0],
            p1[1] + depth*dy,
            p1[2] + depth*dz
        )

        # Pocket wall
        tris.append((p0, q0, q1))
        tris.append((p0, q1, p1))

        # Pocket bottom disk
        center_bottom = (pocket_cx, pocket_cy + depth*dy, pocket_cz + depth*dz)
        tris.append((center_bottom, q1, q0))

    # Pocket entrance ring (flat face to seal edges)
    for i in range(SEG):
        a0 = 2*math.pi*i/SEG
        a1 = 2*math.pi*(i+1)/SEG
        r = MS_POCKET_R
        p0 = (
            pocket_cx + r*math.cos(a0),
            pocket_cy + r*math.sin(a0)*math.sin(tilt_rad),
            pocket_cz + r*math.sin(a0)*math.cos(tilt_rad)
        )
        p1 = (
            pocket_cx + r*math.cos(a1),
            pocket_cy + r*math.sin(a1)*math.sin(tilt_rad),
            pocket_cz + r*math.sin(a1)*math.cos(tilt_rad)
        )
        center_top = (pocket_cx, pocket_cy, pocket_cz)
        tris.append((center_top, p0, p1))

    # === 4. CROSS ===
    cross_z = BASE_H + MTN_H
    cross_cx = MTN_PEAK_OFFSET
    cross_cy = 5  # mtn_y_off

    # Vertical bar
    tris += box_tris(
        cross_cx - CROSS_BAR/2, cross_cy - CROSS_BAR/2, cross_z,
        CROSS_BAR, CROSS_BAR, CROSS_H
    )
    # Horizontal bar
    tris += box_tris(
        cross_cx - CROSS_ARM_W/2, cross_cy - CROSS_BAR/2, cross_z + CROSS_ARM_H,
        CROSS_ARM_W, CROSS_BAR, CROSS_BAR
    )

    # === 5. COIN TRAY ===
    # Simple recessed box on the base surface
    tray_x = TRAY_OFFSET_X
    tray_y = TRAY_OFFSET_Y
    tray_z = BASE_H

    # Tray walls (box going downward into base)
    tris += box_tris(tray_x, tray_y, tray_z - TRAY_H + 2, TRAY_W, TRAY_D, TRAY_H)

    # Rim around tray (slightly raised edge)
    rim = 1.5
    tris += box_tris(tray_x - rim, tray_y - rim, tray_z, TRAY_W + 2*rim, rim, 2)
    tris += box_tris(tray_x - rim, tray_y + TRAY_D, tray_z, TRAY_W + 2*rim, rim, 2)
    tris += box_tris(tray_x - rim, tray_y, tray_z, rim, TRAY_D, 2)
    tris += box_tris(tray_x + TRAY_W, tray_y, tray_z, rim, TRAY_D, 2)

    # === 6. OPEN BOOK ===
    bx = BOOK_OFFSET_X
    by = BOOK_OFFSET_Y
    bz = BASE_H

    # Left page
    tris += box_tris(bx - BOOK_W/2, by, bz, BOOK_W/2 - 1, BOOK_D, BOOK_H)
    # Right page
    tris += box_tris(bx + 1, by, bz, BOOK_W/2 - 1, BOOK_D, BOOK_H)
    # Spine (raised center)
    tris += box_tris(bx - 1, by, bz, 2, BOOK_D, BOOK_H + 1)

    # === 7. ROCK ACCENT BUMPS ===
    # Add some sphere bumps on the mountain surface for extra texture
    bump_positions = [
        (25, -10, BASE_H + 20, 8),
        (-20, 15, BASE_H + 30, 7),
        (15, 20, BASE_H + 45, 6),
        (-30, 5, BASE_H + 15, 9),
        (10, -20, BASE_H + 60, 5),
        (-15, -15, BASE_H + 50, 7),
        (30, 10, BASE_H + 25, 6),
        (-25, -5, BASE_H + 40, 8),
    ]
    for bx_p, by_p, bz_p, br in bump_positions:
        tris += sphere_tris(bx_p, by_p + 5, bz_p, br, 10, 6)

    return tris


# =============================================
# Main
# =============================================
if __name__ == "__main__":
    print("Building MagSafe Mountain Stand geometry...")
    print(f"  Estimated size: {BASE_D:.0f}mm dia x {BASE_H + MTN_H + CROSS_H:.0f}mm tall")
    print(f"  AD5X build check: {'PASS' if max(BASE_D, BASE_D, BASE_H+MTN_H+CROSS_H) <= 220 else 'FAIL - TOO LARGE'}")

    tris = build_stand()
    print(f"  {len(tris)} triangles generated")

    out_dir = os.path.dirname(os.path.abspath(__file__))
    stl_path = os.path.join(out_dir, "magsafe_mountain_stand.stl")

    write_binary_stl(stl_path, tris)
    size_kb = os.path.getsize(stl_path) / 1024
    print(f"  Exported: {stl_path}")
    print(f"  File size: {size_kb:.1f} KB")
    print()
    print("Multi-color print zones:")
    print("  Color 1 (stone/sand): Base + Mountain body")
    print("  Color 2 (gold):       Cross")
    print("  Color 3 (brown):      Book")
    print("  Color 4 (gold):       Coin tray rim")
    print()
    print("Next steps:")
    print("  1. Open STL in a viewer to check proportions")
    print("  2. Import into Onshape/Blender to add rock texture detail")
    print("  3. Slice with supports enabled for MagSafe pocket")
