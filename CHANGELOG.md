# Changelog

Changes a caller has to act on. Full notes for every version are on the
[releases page](https://github.com/JuliaSMLM/MicroscopePSFs.jl/releases).

## [0.5.7] - 2026-11-02

### Changed

- `PupilFunction(zc)` reads `ZernikeCoefficients` magnitudes as absolute. The amplitude is
  `Σⱼ mag[j] Zⱼ`; it used to be `1 + Σⱼ mag[j] Zⱼ`. A caller holding offsets from unit
  amplitude must add 1 to `mag[1]`, and `VectorPSF(...; base_zernike)` and
  `ScalarPSF(...; zernike_coeffs)` read the same way.
  - `mag[5] = 0.06` on an all-zero vector now gives `0.06 Z₅`, not `1 + 0.06 Z₅`. It raises
    no error, and the normalized pupil is wrong.
  - All-zero magnitudes now give a zero pupil. `PupilFunction(zc)` returns it with no
    error; `ScalarPSF(...; zernike_coeffs)` and `VectorPSF(...; base_zernike)` refuse it.
    A phase-only pupil needs `mag[1] = 1`.
  - Coefficients from `ZernikeCoefficients(n)` already carry `mag[1] = 1` and need no edit
    for phase aberrations. Magnitude terms set on them now count at full size against the
    piston, where the added 1 used to halve them. To keep their earlier size relative to
    the piston, halve them.
