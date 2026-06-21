# MolSurfComp → Python porting contract (READ FIRST)

Goal: faithful Python port of the MATLAB package at `/home/ishitani/ext/MolSurfComp`
that reads an xyzr file and writes a triangulated PLY of the Solvent Excluded
Surface (SES). You are porting ONE module. Follow this contract exactly so all
modules integrate. The pipeline (in `MolSurfComp.m` / `Core/data_MolSurf.m`):

```
read_xyzr -> interstructure -> data_I_Cir -> data_Seg_Pat -> data_ext
          -> data_SESsphpat (convex) + SESconcavepat (concave) + data_SEStorpat (toroidal)
          -> accumulate mesh -> PLY
```

## Indexing convention (CRITICAL — already used by existing code)
To mirror the index-interdependent MATLAB source we use **1-based indexing**:
- **Entity arrays** carry a dummy slot at index 0; real entities are 1..N.
- **Coordinate arrays** (centers, `I` points): dummy ROW 0, real columns 0,1,2 = x,y,z.
  Access whole point as `I[s]` (a length-3 vector); `C[i]` is atom i's center.
- **Record matrices** (circle, segment, ncrasegment, segment0, direction, Iijk,
  loops, patches, satom, Ii, circleindex, ...): dummy ROW 0 **and** dummy COLUMN 0,
  so MATLAB `X(r, a:b)` → Python `X[r, a:b+1]` (same literal start; +1 on slice end
  because Python slice end is exclusive). MATLAB `X(r,c)` → `X[r, c]`.
- Convert to 0-based ONLY for final PLY face indices (handled by MeshState).
- MATLAB `mod(m, N) + 1` (1-based circular next) when `m` is 1-based stays
  `m % N + 1` in this scheme (do NOT shift). Be careful porting `mod`.

## Existing files you MUST read and match (already implemented, do not change)
- `molsurf/geom.py`   — `Geom(M, centers[(M+1,3)], R[(M+1,)])`, `read_xyzr`.
- `molsurf/params.py` — `Para(radius_probe, mesh_size)`.
- `molsurf/intersection.py` — `Inter(M_int[(M+1,kmax+1)], num_int[(M+1,)])`,
  `interstructure(geom, Rp)`. `M_int[i, 1:num_int[i]+1]` = sorted neighbor atoms.
- `molsurf/mathutil.py` — `norm, unit, cross, triple, arc_angle(u,v,n,direct=1),
  orthogonalvectors(n)->(v1,v2), circlecenter(c1,c2,r1,r2)`. REUSE these; the
  MATLAB local `alpha(direct,u,v,n)` == `arc_angle(u,v,n,direct)`; the concave
  `alpha(u,v,n)` (no direct) == `arc_angle(u,v,n,1)`.
- `molsurf/sas.py` — `DataI`, `DataCir`, `data_I_Cir(geom, inter, Rp)` (DONE).
  Field layouts (1-based):
  - `DataI.I[(nI+1,3)]`, `.nI`, `.direction[(nI+1,4)]` cols 1..3,
    `.Iijk[(nI+1,4)]` cols 1..3 = i,j,k, `.Ii[(M+1,*)]`, `.In[(M+1,)]`,
    `.I_circle[(M+1,kmax+1,*)]`, `.I_circle_num[(M+1,kmax+1)]`,
    `.high_I[(nI+1,)]` (1=non-singular triangle, 0=singularity candidate),
    `.hightvalue[(nI+1,)]`.
  - `DataCir.circle[(ncircle+1,10)]`: cols `[_, i, j, A(3)=3:6, n(3)=6:9, r=9]`,
    `.ncircle`, `.circleindex[(M+1,11)]`, `.ncircleindex[(M+1,)]`.
- `molsurf/mesh_state.py` — `MeshState` accumulator with
  `add_patch(P, T, face_normals)` where `P[1:]` are real (3,) points (dummy P[0]),
  `T` = iterable of 1-based (a,b,c) index triples into P, `face_normals` aligned
  with T (or None). USE THIS as the single mesh sink for ALL patch families.

## Data structures produced by data_Seg_Pat (agent A defines; others consume)
- `DataSeg.segment[(nseg+1,6)]` cols `[_, i, j, p1, p2, direct]` (p1,p2 index into I; direct ∈ {1,-1}).
- `DataSeg.ncrasegment[(nseg+1,9)]` cols `[_, n(3)=1:4, A(3)=4:7, r=7, radian=8]`.
- `DataSeg.satom[(M+1,*)]`, `DataSeg.nsatom[(M+1,)]`.
- `DataLoop.loops[(nloops+1,smax+1)]` (entries = global segment indices),
  `DataLoop.loopsize[(nloops+1,)]`, `DataLoop.nloops`, `DataLoop.loops_index[(M+1,3)]` cols [_,start,end].
- `DataPat.patches[(npatch+1,pmax+1)]` (entries signed: +k loop index local to atom,
  -k circle index local to atom), `.patchesize[(npatch+1,)]`, `.npatches`,
  `.patches_index[(M+1,3)]` cols [_,start,end], `.patch_atom[(npatch+1,)]`.
- `DataAV` mutable dataclass with floats Acsas,Acses,Aesas,Aeses,Vcsas,Vcses,Vesas,Veses
  (+ running V_cSES,V_eSES,V_cSAS,V_eSAS for concave/toroidal subtraction). Pass by reference.

## mesh_sphpat interface (agent C defines; convex & concave drivers call it)
```
mesh_sphpat(state, c_sphere, r_sphere, loops, loopsize, segment0, circle0,
            patches, patchesize, *, Rp, d, Rj=None)
```
- `state`: MeshState. `c_sphere`: np(3). `r_sphere`: float (= R_i+Rp for convex, = Rp for concave).
- `loops[k][1:loopsize[k]+1]` = indices into `segment0` (the segments of local loop k).
- `segment0[(*,12)]` cols `[_, c(3)=1:4, n(3)=4:7, r=7, spoint(3)=8:11, angle=11]`,
  n points OUTWARD from the sphere, angle is the arc radian.
- `circle0[(*,8)]` cols `[_, c(3)=1:4, n(3)=4:7, r=7]` (or None / shape (1,8) dummy if no circles).
- `patches[1:patchesize+1]`: signed; `+k` → use `loops[k]`, `-k` → use `circle0[k]`.
- `Rj`: optional dict/array mapping global segment index → neighbour atom VdW radius
  (only used by the convex SAS→SES near-cusp arc refinement; pass None for concave → rj=Rp).
- Internally: `tolerance = 0.8*min(d, r_sphere)`; for convex (`r_sphere>Rp`) the
  boundary and sampled points are contracted by factor `(r_sphere-Rp)/r_sphere` toward
  c_sphere (this is the SAS→SES contraction in loop_division/circle_division/map_sphere).
  For concave (`r_sphere==Rp`) no contraction.
- Normals via compute_NV: convex `arg_NV=+1` (outward = away from atom center),
  concave `arg_NV=-1` (outward = toward probe center). Then `state.add_patch(P,T,NV)`.

## Conventions for all agents
- numpy only. Python 3.11. Reuse `molsurf/mathutil.py` helpers; clip all `acos`
  arguments to [-1,1]; use `np.floor(...).astype(int)` for cell indices (never `int()`
  on negatives); guard `sqrt(max(x,0))` where MATLAB silently took real part.
- Recursion: the advancing front recurses; call `sys.setrecursionlimit(100000)` at
  import in meshing.py.
- Write clean, faithful, line-by-line ports. Keep MATLAB function names as Python
  function names where reasonable. Add a short comment citing the MATLAB source file.
- Do NOT port visualization (`visu*`, figures) or the analytic text writers
  (`output_SES_*patches`, `output_SAS_patches`) — those are out of scope. Port only
  the geometry + mesh-accumulation paths.
- Each agent: after writing, run a smoke import (`.venv/bin/python -c "import molsurf.<mod>"`)
  and any quick self-test you can, and report results.
