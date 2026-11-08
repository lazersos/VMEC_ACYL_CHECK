# Compare GetBcyl with curl(GetAcyl)

Small Fortran program linked to an existing LIBSTELL build. It reads a WOUT
file from the command line, samples physical cylindrical vector components on
a uniform (R,phi,Z) grid, and compares GetBcyl with a centered finite-difference
curl of GetAcyl. No rescaling or sign correction is applied to either result.

## Build

Use the same Fortran compiler/version and compatible flags as the LIBSTELL
build. The existing modules must expose the GetBcyl and GetAcyl generic
interfaces (including their WOUT implementations).

Set the directory containing stel_kinds.mod, read_wout_mod.mod, and
vmec_utils.mod, the archive path, and the archive's link dependencies:

```sh
make FC=gfortran \
  MODDIR=/absolute/path/to/fortran/modules \
  LIBSTELL=/absolute/path/to/libstell.a \
  LDLIBS="$(nf-config --flibs) -llapack -lblas"
```

The dependency example is illustrative: reuse the libraries and link order
from a working STELLOPT executable. A particular build may also need HDF5,
FFTW, OpenMP flags, or an MPI compiler wrapper such as mpifort. Modules and
archives may live in different directories; MODDIR and LIBSTELL are independent.
For multiple module directories, add `CPPFLAGS='-I/another/module/directory'`.
Paths used in Makefile variables should not contain spaces.

You can put local settings in config.mk instead:

```make
FC = gfortran
MODDIR = /absolute/path/to/modules
LIBSTELL = /absolute/path/to/libstell.a
LDLIBS = -L/path/to/netcdf/lib -lnetcdff -lnetcdf -llapack -lblas
```

LIBSTELL itself must have been compiled with -g for library source line
information, and with -fcheck=all for checks inside its routines. On macOS:

```sh
make dsym
```

This runs dsymutil on the final executable while the library is available.

## Run

```sh
./compare_wout /path/to/wout_example.nc
./compare_wout /path/to/wout_example.nc grid.nml
```

Defaults are 33 x 33 x 33 samples, a rectangular box estimated from sampled
boundary Fourier series, and one toroidal field period. The first and last
phi planes are duplicates. Only nphi-1 unique planes are evaluated and counted.
Use full_torus=.true. to cover 2*pi instead. Angles are in radians and lengths
in meters. Standard output prints progress and comparison statistics.

grid.nml controls dimensions, bounds, phi0, full_torus, the accepted s range,
and the output filename. Its commented R/Z limits are examples; leave them
commented to use the WOUT-derived defaults. All arrays are allocated by the
program, and the curl routine is in a module with an explicit interface.

Only samples with successful return codes, finite values, and
`s_min < s <= s_max` are accepted. B is sampled first to seed GetAcyl's inverse
map. GetAcyl is evaluated only where that inverse map succeeded. A comparison
requires valid A at the center and all six stencil neighbors, plus a valid
GetBcyl result at the center. No one-sided differences or exterior zero-filled
values are used. Radial and vertical grid edges are omitted. The exact magnetic
axis is excluded even when s_min=0.

## Output and interpretation

comparison.csv contains only valid comparison points, with columns:

```
R,phi,Z,s,Br_getBcyl,Bphi_getBcyl,Bz_getBcyl,
Br_curlA,Bphi_curlA,Bz_curlA,dBr,dBphi,dBz,abs_error,relative_error
```

Differences are curl(A)-B. abs_error is the vector norm in tesla;
relative_error is that norm divided by the local reference field norm.
Duplicate endpoint rows are not written. The program exits unsuccessfully
if there are no valid comparisons or the reference norm vanishes.

The relative RMS is sqrt(sum(|curl(A)-B|^2)/sum(|B|^2)). The RMS vector error
is sqrt(sum(|curl(A)-B|^2)/N). These are unweighted grid-point statistics, not
cylindrical volume integrals. The fitted scale is
sum(curl(A).B)/sum(|B|^2): ideally near 1; a value near 2*pi or -1 suggests a
normalization or sign error. It is reported only and never used to alter data.

Repeat with doubled interval counts, e.g. nr=nphi=nz=65 after a 33-point run.
The centered stencil is second order for smooth fields. GetBcyl and GetAcyl
may use different radial interpolants, so agreement can reach an interpolation
floor; their finite-resolution representations need not agree exactly.
The compared point set also changes on refinement, so the aggregate RMS ratio
alone is not a rigorous convergence measurement on a fixed physical domain.

## GetAcyl_WOUT prerequisites

This project tests the installed GetAcyl; it does not replace or patch it.
In the alternative gauge discussed, WOUT fluxes give covariant components
As=-isigng*phipf*lambda/(2*pi), Au=isigng*phi/(2*pi), and
Aphi_cov=-chi/(2*pi). Unlike phi and phipf, the WOUT chi profile already
carries VMEC's internal orientation: eqfor.f integrates the internal phip
times iotas and writes chi=2*pi*chi1. Do not multiply chi by isigng again.
With isigng=-1, that extra factor would reverse the axisymmetric poloidal
field while leaving the toroidal field unchanged.
The physical toroidal component is Aphi_cov/R
only after accounting for the Rphi and Zphi contributions in the transformation.

Check your implementation for the 1/(2*pi) normalization and for sflx/uflx
being INTENT(INOUT) if they are read as initial guesses. Some versions of
tosuvspaceLambda write all three optional outputs unconditionally at entry;
in those versions GetAcyl must pass LAM, LU, and LV, even if only LAM is used.
Older text WOUT readers may not allocate chi; use a NetCDF WOUT with the
required profiles, or adapt GetAcyl to reconstruct its poloidal flux.
The comparison intentionally leaves these behaviors visible.

## Plot a toroidal plane

The included plot_comparison.py requires NumPy and Matplotlib:

```sh
python plot_comparison.py comparison.csv --phi-index 8 -o phi8.png
python plot_comparison.py comparison.csv --list-phi
```

The default figure has rows for BR, Bphi, BZ, and field magnitude, and columns
for GetBcyl, curl(GetAcyl), and relative error. Indices are one based among
sorted planes present in the CSV. Component errors are normalized by the
reference vector magnitude; the magnitude row shows vector relative error.
Missing cells remain blank. Use --component magnitude for a three-panel figure.

## Verification

```sh
make check MODDIR=/path/to/LIBSTELL/modules
```

The analytic curl test needs only stel_kinds.mod, not libstell.a. It exercises
all three components on a nonaxisymmetric periodic field, verifies second-order
toroidal convergence, checks duplicated endpoints and invalid stencils, and
runs with bounds checking under the default GNU flags.

See VALIDATION.txt for checks performed when this package was created.
