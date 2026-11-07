# DUO Aquarium — 3-bit demonstration (v0.1)

This is an **experimental teaching and research toy**, not a physical simulation or a validated thermodynamic model. It is a **standalone C program** that does not require GMP or any original Dynamic Unary source code.

## Build on Fedora

For the headless version, GCC and make suffice:

```bash
make aquarium
./aquarium --headless --steps 100000 --seed 42
```

For the windowed version, install the X11 development package if needed:

```bash
sudo dnf install gcc make libX11-devel
make aquarium-gui
./aquarium-gui
```

If running through SSH without graphical forwarding, choose `--headless`. GUI controls: **Space** pause/resume; **N** one step while paused; **Q/Esc** quit.

Both can be compiled using `make` when X11 development headers are present. Binaries are built on **your own computer** for compatibility.

## Three-bit cycles

```
A: 000 -> 110 -> 011 -> 010 -> 000
B: 111 -> 001 -> 100 -> 101 -> 111
```

There are exactly 8 states and 2 precomputed orbits. No DUE transform is computed at runtime. Each object stores an orbit index, an integer phase (0..3), a fractional clock, a spin rate, position and spatial velocity.

## Rules

1. All objects rotate **forward** through their own current orbit; spin rate controls how many phase transitions occur per second.
2. Objects move in a 2D rectangular arena and bounce off walls. Use `--periodic` to wrap around opposite edges instead.
3. When overlapping objects are moving toward each other, they have an equal-mass elastic spatial collision.
4. Compute `C = state_A XOR state_B`; **both objects enter the orbit and phase containing C**, retaining their *independent* fractional clocks and their own spin rates.
5. Exchange **1/8 of the signed difference** in spin rate: `transfer = (spin_A - spin_B)/8`, `spin_A -= transfer`, `spin_B += transfer`. This conserves the **sum of spin rates**, which we use as a heat-like *proxy*, not real physical energy. They do not instantly average.
6. Collisions can discard prior binary information; the overall collision rule is **not reversible**. Native DU cycles remain reversible individually.

## Data

By default, files are written to the current directory:

- `aquarium_events.csv`: one row per collision, with tick, object IDs, binary states before the collision, XOR result, new orbit/phase, spin rates before/after, and impact location. Values 0..7 are **decimal encodings of 3-bit patterns**.
- `aquarium_samples.csv`: one row every 60 simulation ticks, with population in each orbit, total/mean/spread/min/max spin rate, and counts of all eight binary patterns. An initial snapshot at tick zero is also included.

Use `--prefix experiment_01` to select the filenames; use `--no-log` to disable files. To repeat a run, use the same seed, object count, boundary mode, and steps. Use separate prefixes to retain previous logs.

## Examples

```bash
./aquarium --headless --seed 1234 --steps 60000 --prefix test1
./aquarium --headless --seed 1234 --steps 60000 --prefix test1_repeat
cmp test1_samples.csv test1_repeat_samples.csv
./aquarium-gui --periodic --seed 98 --objects 24 --prefix periodic98
```

The simulation uses **60 discrete time steps per simulated second**. The X11 viewer runs close to real time; headless mode runs as fast as the CPU permits. GUI window size is 800 × 520 plus a status area. Rings around objects show a two-band visual indicator of their spin rate, not literal temperature.

## Questions worth studying

- Do orbit-A and orbit-B populations approach a regular distribution?
- Do spin rates converge? Our chosen heat exchange rule directly pushes them toward one another, so such convergence is *expected*, not emergent thermodynamics.
- Do distinct local phases persist even when every object belongs to the same orbit?
- Does changing the random seed or arena boundary conditions alter long-run statistics?
- How do outcomes differ when collisions only exchange phase versus forcing the same resulting phase?

## Limits of this model

The present model does not equate rotation rate with thermodynamic energy, include physical temperature, enforce microscopic reversibility of collisions, guarantee exact event chronology independent of the 60 Hz step size, or implement continuous collision detection. The simulation has no claims of emergent physics. It is intended to let us explore **informational interactions** with clearly specified and reproducible rules.
