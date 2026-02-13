---
name: pic-to-part
description: Analyze a reference photo or image of a physical object and generate a 3D printable model from it. Works with photos, screenshots, sketches, or any visual reference.
argument-hint: [key measurements, e.g. "40mm wide, 15mm tall"]
allowed-tools:
  - Write
  - Edit
  - Read
  - Bash
---

# Picture to Part Generator

You are analyzing a reference image provided by the user and generating a 3D printable model from it.

## Workflow

1. **Analyze the image**: Identify the overall shape, features, symmetry, holes, curves, and any visible details.
2. **Ask for scale** (if not provided): You need at least ONE real-world measurement to set the scale. Ask: "What's one key dimension I can use for scale? (e.g., overall width = 40mm)"
3. **Estimate proportions**: Use the reference dimension to calculate all other dimensions from the image proportions.
4. **Choose approach**:
   - Simple geometric shapes → OpenSCAD
   - Complex curves/fillets → CadQuery (Python)
5. **Generate the model** in `parts/<part-name>/`
6. **Run the Python STL exporter** to produce the .stl file immediately.
7. **Summarize** what you modeled and list any assumptions you made about dimensions.

## Analysis Checklist

When looking at the image, identify:
- [ ] Overall bounding box (L x W x H proportions)
- [ ] Symmetry (bilateral, radial, none)
- [ ] Primitive shapes (boxes, cylinders, spheres, cones)
- [ ] Boolean operations (holes, cutouts, slots)
- [ ] Fillets and chamfers (rounded edges)
- [ ] Repeating patterns (arrays, grids)
- [ ] Thin features (walls, ribs, flanges)
- [ ] Interface features (snap fits, screw holes, slots)

## Proportion Estimation

Given one known dimension D_known on axis A:
1. Measure the pixel dimensions of D_known in the image
2. Calculate scale: mm_per_pixel = D_known / pixel_count
3. Apply to all other visible features
4. Round to reasonable precision (0.1mm for small parts, 0.5mm for larger)

## Output Requirements

- Create `parts/<part-name>/` directory
- Generate `measurements.md` with estimated dimensions and assumptions
- Generate the .scad or .py model file
- Generate `export_stl.py` if using OpenSCAD approach
- Run the export to produce the .stl file
- Tell the user:
  - What assumptions were made
  - Which dimensions are estimated vs. provided
  - Suggest which measurements to verify with calipers

## Printer: Flashforge Adventurer 5X (AD5X)

All generated parts are validated against these specs:
- **Build volume**: 220 x 220 x 220mm — WARN if part exceeds this
- **Nozzle**: 0.4mm — min wall = 0.8mm, min feature = 0.4mm
- **Layer heights**: 0.12mm (fine detail), 0.2mm (normal), 0.28mm (draft)
- **Materials**: PLA, PETG, TPU 95A, PLA-CF, PETG-CF
- **Multi-color**: 4 colors via IFS — note color zones in the model if applicable
- **Hole compensation**: +0.2mm to nominal diameter for 0.4mm nozzle

Always include recommended print settings (layer height, infill, supports, orientation) based on the part geometry.

## Important Notes

- Always be transparent about estimated vs. measured dimensions
- Err on the side of slightly larger tolerances for estimated parts
- If the image is ambiguous, ask the user before guessing
- For parts that mate with other components, ask for critical interface dimensions
- Suggest the user verify key dimensions before printing
