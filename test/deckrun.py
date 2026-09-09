"""Run a shipped example deck and read its output, for verification tests.

Both solver-level verification tests - the RMF antenna in vacuum and the RMF
boundary condition in vacuum - do the same three things: copy a deck and its
mesh into a scratch directory, preprocess and run it, and read the .vtu frames
back in time order. That is all this module is.

Two details in here are not obvious and were each got wrong once:

  * frames must be ordered by their INDEX, not by filename. A plain sort gives
    _0, _1, _10, _11, _12, _2, ... which silently scrambles time order, and any
    test that takes "the later frames" or measures a rotation between
    consecutive frames then reads a shuffled sequence.

  * the solver does not exit non-zero on every kind of bad solution, so the log
    has to be searched for the diagnostics it prints instead.
"""

import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, 'scripts')

try:
    import numpy as np
except ImportError:  # pragma: no cover - callers skip themselves without numpy
    np = None

try:
    import vtu
except ImportError:  # pragma: no cover
    vtu = None

# Diagnostics the solver prints when the solution has gone bad. It does not
# always exit non-zero on these.
BAD = ('NaN', 'Negative pressure', 'Negative density', 'Negative electron',
       'Negative ion')


def run(example_dir, pin_name, workdir, aux_files=(), binary=None, timeout=3600):
    """Preprocess and run one deck in its own directory; return that directory.

    aux_files are copied alongside the deck - meshes, mostly, which the deck
    names relative to the working directory.
    """
    case = os.path.join(workdir, os.path.splitext(pin_name)[0])
    os.makedirs(case, exist_ok=True)
    for name in (pin_name,) + tuple(aux_files):
        shutil.copy(os.path.join(example_dir, name), case)

    env = dict(os.environ, PYTHONPATH=SCRIPTS)
    pre = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, 'wxinpparse.py'), '-i', pin_name],
        cwd=case, env=env, capture_output=True, text=True, timeout=300)
    if pre.returncode != 0:
        raise RuntimeError(f'preprocessing {pin_name} failed:\n{pre.stdout}{pre.stderr}')

    inp = os.path.splitext(pin_name)[0] + '.inp'
    run_ = subprocess.run([binary or apollo_binary(), '-i', inp],
                          cwd=case, capture_output=True, text=True, timeout=timeout)
    log = run_.stdout + run_.stderr
    if run_.returncode != 0:
        raise RuntimeError(f'solver exited {run_.returncode}:\n{log[-2000:]}')
    for bad in BAD:
        if bad in log:
            raise RuntimeError(f'solver reported "{bad}":\n{log[-2000:]}')
    return case


def apollo_binary():
    """The solver to run, absolute so that a per-case cwd cannot break it."""
    if os.environ.get('APOLLO_BIN'):
        return os.path.abspath(os.environ['APOLLO_BIN'])
    for variant in ('build-opt', 'build-debug'):
        path = os.path.join(REPO, 'src', variant, 'apollo')
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    return None


def frame_index(name):
    return int(name.rsplit('_', 1)[1].split('.')[0])


def frames(case, minimum=4):
    """[(index, data)] for every .vtu in `case`, ordered by frame index."""
    names = sorted((f for f in os.listdir(case) if f.endswith('.vtu')),
                   key=frame_index)
    if len(names) < minimum:
        raise RuntimeError(f'expected at least {minimum} output frames in {case}, '
                           f'got {len(names)}')
    return [(frame_index(n), vtu.read(os.path.join(case, n))) for n in names]


def nodal(data, components, neq=18):
    """Nodal coordinates and the named state components, as flat arrays.

    For P1 on triangles the three nodes are the cell vertices, and the cell-data
    arrays are laid out node-major: solutiondg.(node*neq + component).
    """
    points, cells = data['Position'], data['connectivity'].reshape(-1, 3)
    xs, ys = [], []
    out = {c: [] for c in components}
    for node in range(3):
        vid = cells[:, node]
        xs.append(points[vid, 0])
        ys.append(points[vid, 1])
        for c in components:
            out[c].append(data[f'solutiondg.{node * neq + c}'])
    return (np.concatenate(xs), np.concatenate(ys),
            {c: np.concatenate(v) for c, v in out.items()})
