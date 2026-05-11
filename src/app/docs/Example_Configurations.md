# Example Configurations

All configuration files live in the `configs/` directory and are passed to the executable
via the `-c` flag:

```bash
mtrsim -c configs/<file>.json -o /path/to/output
```

Every file uses `"seed": 42` so results are directly comparable to one another and to the
default run.

---

## Image Size

The output PNG is one pixel per voxel.  Width and height are determined entirely by the
physical volume dimensions and the voxel spacing:

```
image width  (px) = round(xLen / dx)
image height (px) = round(yLen / dy)
```

There are two independent ways to change the image size:

- **Change `dx` / `dy`** — adjusts pixel density without altering the physical volume.
  Halving the spacing doubles each dimension and quadruples pixel count, but also doubles
  the Cholesky matrix size in each direction, so runtime increases significantly.
- **Change `xLen` / `yLen`** — adjusts the physical area simulated.  Doubling the volume
  at the same spacing doubles each pixel dimension and shows more MTR regions, but the
  per-MTR size in pixels stays the same.

| Config file | Image size | vs. default | Description |
|---|---|---|---|
| `default.json` | 1905 × 635 px | — | Baseline — matches `simulate_MTRs.m` defaults |
| `image_low_res.json` | 953 × 318 px | ¼ pixels, ~8× faster | Same physical volume (`dx=dy=0.04 mm`), coarser voxels |
| `image_high_res.json` | 3810 × 1270 px | 4× pixels, ~4–8× slower | Same physical volume (`dx=dy=0.01 mm`), finer voxels |
| `image_large_volume.json` | 3810 × 1270 px | 4× pixels, similar cost/voxel | Double the physical area (`xLen=76.2, yLen=25.4 mm`), MTR physical sizes unchanged |

### Key values

**`image_low_res.json`** — halved resolution, fast preview:
```json
"dx": 0.04,  "dy": 0.04,
"xLen": 38.1, "yLen": 12.7
```

**`image_high_res.json`** — doubled resolution, higher detail:
```json
"dx": 0.01,  "dy": 0.01,
"xLen": 38.1, "yLen": 12.7
```

**`image_large_volume.json`** — larger physical scan area, same voxel spacing:
```json
"dx": 0.02,  "dy": 0.02,
"xLen": 76.2, "yLen": 25.4
```

---

## MTR Region Size and Shape

MTR size and morphology are controlled by `thetaList`.  This is a 2D array with one row
per latent Gaussian field (`numComponents − 1` rows, always 2 for the default
3-component system) and three columns for the x, y, z spatial directions.

The underlying correlation function is exponential:

```
rho(distance, theta) = exp( −distance / theta )
```

Correlation drops to ~37% at a lag equal to `theta`.  All `theta` values are in mm,
the same units as `dx` and `xLen`.

The default values produce the needle-like MTR morphology characteristic of titanium:
large `theta_y` (long in the rolling/elongation direction) and small `theta_x`
(narrow across it).

All MTR size configs below use the default image size (1905 × 635 px).

| Config file | theta_x (mm) | theta_y (mm) | Morphology |
|---|---|---|---|
| `default.json` | 0.10 / 0.08 | 0.45 / 0.37 | Baseline elongated needle-like MTRs |
| `mtrs_large.json` | 0.20 / 0.16 | 0.90 / 0.74 | ~2× physically larger MTRs |
| `mtrs_small.json` | 0.04 / 0.03 | 0.18 / 0.15 | Fine-grained, many small MTRs |
| `mtrs_isotropic.json` | 0.30 / 0.25 | 0.30 / 0.25 | Blocky, equiaxed — no preferred elongation direction |
| `mtrs_very_elongated.json` | 0.05 / 0.04 | 1.50 / 1.20 | Long thin streaks, ~30:1 aspect ratio |

*Each row in the theta columns above shows Gaussian Y₁ / Gaussian Y₂.*

### Key values

**`mtrs_large.json`** — doubled theta, physically larger regions:
```json
"thetaList": [
  [0.20, 0.90, 0.20],
  [0.16, 0.74, 0.16]
]
```

**`mtrs_small.json`** — reduced theta, fine-grained texture:
```json
"thetaList": [
  [0.04, 0.18, 0.04],
  [0.03, 0.15, 0.03]
]
```

**`mtrs_isotropic.json`** — equal theta in all directions, blocky morphology:
```json
"thetaList": [
  [0.30, 0.30, 0.30],
  [0.25, 0.25, 0.25]
]
```

**`mtrs_very_elongated.json`** — exaggerated needle-like morphology:
```json
"thetaList": [
  [0.05, 1.50, 0.05],
  [0.04, 1.20, 0.04]
]
```

---

## Combining Image Size and MTR Size

The two controls are fully independent.  To produce, for example, a high-resolution
image of large MTRs, combine the relevant values from each group in a single config file:

```json
{
  "dx": 0.01,
  "dy": 0.01,
  "xLen": 38.1,
  "yLen": 12.7,
  "thetaList": [
    [0.20, 0.90, 0.20],
    [0.16, 0.74, 0.16]
  ],
  "odfInputPath": "data/simulation_ODF.h5",
  "seed": 42
}
```

---

## Runtime Guidance

| Scenario | Approximate relative cost |
|---|---|
| Default (`dx=0.02`, `xLen=38.1`) | 1× |
| Low-res (`dx=0.04`, same volume) | ~0.1× |
| High-res (`dx=0.01`, same volume) | ~8–16× |
| Large volume (`xLen=76.2`, `dx=0.02`) | ~8–16× |
| MTR size changes only | Same as default |

Runtime is dominated by the Cholesky decompositions, which scale as O(nx²) and O(ny²).
Halving `dx` doubles `nx` and roughly quadruples the Cholesky cost in the x direction.
Changes to `thetaList` alone do not affect runtime.
