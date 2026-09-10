#!/usr/bin/env python3
"""Reduce a run's .vtu output to a small, comparable fingerprint.

WHY THIS EXISTS. Apollo's answer depends on things outside the deck: the PETSc
release it was built against (`docs/known-issues.md` §17 is one change that
matters), the compiler, the BLAS. There is no way to notice that from a single
run - the numbers look plausible either way - and the only honest check is to
run the same deck on two builds and compare.

Comparing directly is impractical: a Phase 3 run is up to 1.3 GB, and the two
builds are usually on different machines. This writes a few kilobytes per run
instead, which can be mailed or pasted, and compares two of them.

    # on each machine
    python3 scripts/run_fingerprint.py <results-dir> -o mine.json

    # anywhere
    python3 scripts/run_fingerprint.py --compare theirs.json mine.json

WHAT IT COMPARES, AND WHY THAT SHAPE. Every value of every array, reduced by
statistics that do not depend on cell ordering: count, min, max, and the sum of
the sorted values. Sorting first makes the sum independent of partitioning and
of summation order, which is what lets a one-rank and a two-rank run be compared
at all (`docs/known-issues.md` §15 uses the same device). Two runs that agree
here agree on the multiset of numbers in the file; they could in principle
differ by a permutation, which for this purpose is exactly what we want to
ignore.

WHAT IT WILL NOT DO. It will not tell you a run is "fine". It reports the
largest relative difference it found and where; deciding whether 1e-13 is
round-off and 1e-3 is a bug is your job, and the thresholds it prints are
guides, not verdicts. Where it cannot compare - different frame counts,
different arrays, different sizes - it refuses and says which, rather than
comparing the subset that happens to line up and reporting a small number.
"""

import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, 'test'))

# Arrays that describe the mesh rather than the solution, reported separately so
# a real difference in the fields is not buried among them.
#
# THESE DO NOT ALL BEHAVE THE SAME ACROSS RANK COUNTS, and getting that wrong
# makes the tool useless for the one comparison known-issues 15 is about.
# PETSc writes one <Piece> per rank, and a vertex on a partition boundary
# appears in every piece that touches it - so the concatenated point list is
# LONGER on more ranks (measured: 18303 points on 1 rank, 18588 on 2, same
# mesh) while the cell data is not duplicated and its length is identical
# (11989 in both). A first draft of this file treated any length difference
# here as "the mesh differs" and refused, which reported a partitioning
# difference as two different problems.
TOPOLOGY_ARRAYS = ('Position', 'connectivity', 'offsets', 'types')
# Deliberately its own category: 'Rank' is the rank that owns each cell, so it
# differs between a 1-rank and a 2-rank run by construction. That is not a
# finding.
RANK_ARRAYS = ('Rank',)
MESH_ARRAYS = TOPOLOGY_ARRAYS + RANK_ARRAYS


def _stats(values):
    """Order-independent summary of one array."""
    import numpy as np
    flat = np.asarray(values).ravel()
    if flat.size == 0:
        return {'n': 0, 'min': None, 'max': None, 'sorted_sum': 0.0,
                'abs_sorted_sum': 0.0}
    finite = np.isfinite(flat)
    n_bad = int(flat.size - finite.sum())
    ordered = np.sort(flat[finite])
    return {
        'n': int(flat.size),
        'n_nonfinite': n_bad,
        'min': float(ordered[0]) if ordered.size else None,
        'max': float(ordered[-1]) if ordered.size else None,
        # Summed smallest-first, so the value does not depend on the order the
        # cells happened to be written in.
        'sorted_sum': float(ordered.sum()),
        'abs_sorted_sum': float(np.sort(np.abs(ordered)).sum()),
    }


def fingerprint(directory, frames=None):
    """Fingerprint every .vtu in *directory*, or the frames named."""
    import vtu

    files = sorted(glob.glob(os.path.join(directory, '*.vtu')))
    if not files:
        raise SystemExit('no .vtu files under %s' % directory)
    if frames is not None:
        want = set(frames)
        chosen = [f for i, f in enumerate(files) if i in want]
        if len(chosen) != len(want):
            raise SystemExit(
                'asked for frames %s but the directory holds %d'
                % (sorted(want), len(files)))
        files = chosen

    out = {'source': os.path.abspath(directory), 'frames': []}
    for path in files:
        arrays = vtu.read(path)
        out['frames'].append({
            'file': os.path.basename(path),
            'arrays': {k: _stats(v) for k, v in sorted(arrays.items())},
        })

    log = os.path.join(directory, 'solver.log')
    if os.path.exists(log):
        with open(log) as handle:
            for line in handle:
                # The banner carries the PETSc the binary was built against,
                # which is the whole point of the comparison.
                if line.startswith('Apollo ') and 'PETSc' in line:
                    out['banner'] = line.strip()
                    break
    return out


def _relative(a, b):
    """Relative difference, falling back to absolute when both are tiny."""
    if a is None or b is None:
        return None
    scale = max(abs(a), abs(b))
    if scale == 0.0:
        return 0.0
    # Below this the two are both numerically zero and a ratio is meaningless.
    if scale < 1e-300:
        return 0.0
    return abs(a - b) / scale


def compare(left, right, tolerance):
    """Compare two fingerprints. Returns (worst, rows, refusals, notes)."""
    refusals = []
    notes = []
    if len(left['frames']) != len(right['frames']):
        refusals.append('frame count differs: %d vs %d'
                        % (len(left['frames']), len(right['frames'])))
        return None, [], refusals, notes

    rows = []
    worst = None
    for i, (lf, rf) in enumerate(zip(left['frames'], right['frames'])):
        lkeys, rkeys = set(lf['arrays']), set(rf['arrays'])
        if lkeys != rkeys:
            only_l = sorted(lkeys - rkeys)
            only_r = sorted(rkeys - lkeys)
            refusals.append(
                'frame %d: array names differ (only left: %s; only right: %s)'
                % (i, only_l or '-', only_r or '-'))
            continue
        for name in sorted(lkeys):
            ls, rs = lf['arrays'][name], rf['arrays'][name]
            if ls['n'] != rs['n']:
                if name in TOPOLOGY_ARRAYS:
                    # Expected when the two runs used different rank counts;
                    # the cell data below is still comparable, which is the
                    # whole point. Recorded, not refused.
                    notes.append(
                        'frame %d, %s: %d vs %d points - consistent with a '
                        'different partitioning, not a different mesh'
                        % (i, name, ls['n'], rs['n']))
                    continue
                refusals.append('frame %d, %s: length differs (%d vs %d)'
                                % (i, name, ls['n'], rs['n']))
                continue
            if ls.get('n_nonfinite') or rs.get('n_nonfinite'):
                refusals.append(
                    'frame %d, %s: non-finite values present (%s vs %s)'
                    % (i, name, ls.get('n_nonfinite'), rs.get('n_nonfinite')))
            diffs = {k: _relative(ls[k], rs[k])
                     for k in ('min', 'max', 'sorted_sum', 'abs_sorted_sum')}
            got = max((d for d in diffs.values() if d is not None), default=0.0)
            rows.append((i, name, got, diffs))
            if worst is None or got > worst[2]:
                worst = (i, name, got, diffs)

    rows.sort(key=lambda r: -r[2])
    return worst, rows, refusals, notes


def _report(left, right, tolerance, limit):
    # Every path that cannot compare must reach the caller as a non-zero exit,
    # not only as a printed line. An earlier version of this function printed
    # the refusals and then returned 0 whenever the arrays that DID line up
    # agreed - so a renamed array, a length mismatch, a changed mesh and a file
    # full of NaN all reported "the two runs agree to round-off". Found by
    # mutating a fingerprint and checking the exit status, which is the only
    # thing a CI step or a shell script reads.
    blocked = False
    worst, rows, refusals, notes = compare(left, right, tolerance)

    for side, fp in (('left ', left), ('right', right)):
        print('%s: %s' % (side, fp.get('banner', fp.get('source', '?'))))
    print()

    if refusals:
        blocked = True
        print('*** CANNOT COMPARE:')
        for r in refusals:
            print('***   %s' % r)
        print('*** Refusing rather than comparing whatever lines up: a small')
        print('*** number computed over a subset would read as agreement.')
        print()
        if not rows:
            return 2

    fields = [r for r in rows if r[1] not in MESH_ARRAYS]
    mesh = [r for r in rows if r[1] in MESH_ARRAYS]

    # Is this a same-partitioning comparison or a cross-rank one? It decides
    # what the topology arrays are allowed to say. Under a different
    # partitioning NONE of them is comparable: PETSc duplicates boundary
    # vertices per <Piece>, so Position gets longer; connectivity indexes into
    # that renumbered point list; offsets follow connectivity. Only the cell
    # count and the cell-data multisets survive, and those are exactly what
    # known-issues 15 compares.
    repartitioned = bool(notes) or any(
        r[1] in RANK_ARRAYS and r[2] > tolerance for r in mesh)

    if repartitioned:
        print('note: THE TWO RUNS ARE PARTITIONED DIFFERENTLY (different rank')
        print('      counts, or a different partitioner).')
        for n in notes[:2]:
            print('      %s' % n)
        if len(notes) > 2:
            print('      ... and %d more like it' % (len(notes) - 2))
        print('      Position, connectivity and offsets are therefore NOT')
        print('      comparable and are excluded below - a boundary vertex')
        print('      appears once per piece that touches it, and connectivity')
        print('      indexes into that renumbered list. The cell data is not')
        print('      duplicated, so the solution arrays ARE comparable, which')
        print('      is the comparison known-issues 15 is about.')
        print()
        mesh = []
    else:
        mesh = [m for m in mesh if m[1] not in RANK_ARRAYS]

    if mesh:
        mesh_worst = max(m[2] for m in mesh)
        if mesh_worst > tolerance:
            blocked = True
            print('*** THE MESH DIFFERS (worst %.3e, in %s). Everything below'
                  % (mesh_worst, max(mesh, key=lambda m: m[2])[1]))
            print('*** is then comparing two different problems. Stop here.')
            print()

    print('largest relative differences, solution arrays '
          '(sorted-multiset statistics):')
    print('%-6s %-16s %12s' % ('frame', 'array', 'rel.diff'))
    for i, name, got, _ in fields[:limit]:
        print('%-6d %-16s %12.3e' % (i, name, got))
    if len(fields) > limit:
        print('  ... %d more, all smaller' % (len(fields) - limit))
    print()

    if not fields:
        print('no solution arrays compared')
        return 2

    top = fields[0][2]
    print('worst over %d solution arrays: %.3e  (frame %d, %s)'
          % (len(fields), top, fields[0][0], fields[0][1]))
    if blocked:
        print('NOT COMPARABLE - see the refusals above. The number just')
        print('printed covers only the arrays that lined up, so it is not')
        print('evidence of anything.')
        return 2
    if top <= tolerance:
        print('within the %.1e tolerance: the two runs agree to round-off.'
              % tolerance)
        return 0
    print('ABOVE the %.1e tolerance. These are not the same answer.'
          % tolerance)
    print('This is the signal the comparison exists to produce - do not')
    print('explain it away without finding the cause.')
    return 1


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__.split('\n\n')[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('paths', nargs='+',
                   help='a results directory, or two .json fingerprints '
                        'with --compare')
    p.add_argument('-o', '--output', help='write the fingerprint here')
    p.add_argument('--compare', action='store_true',
                   help='compare two fingerprints written earlier')
    p.add_argument('--frames', type=int, nargs='+',
                   help='fingerprint only these frame indices')
    p.add_argument('--tolerance', type=float, default=1e-10,
                   help='relative difference treated as round-off '
                        '(default 1e-10)')
    p.add_argument('--limit', type=int, default=15,
                   help='how many arrays to list (default 15)')
    args = p.parse_args(argv)

    if args.compare:
        if len(args.paths) != 2:
            p.error('--compare takes exactly two fingerprint files')
        with open(args.paths[0]) as h:
            left = json.load(h)
        with open(args.paths[1]) as h:
            right = json.load(h)
        return _report(left, right, args.tolerance, args.limit)

    if len(args.paths) == 2 and all(os.path.isdir(x) for x in args.paths):
        left = fingerprint(args.paths[0], args.frames)
        right = fingerprint(args.paths[1], args.frames)
        return _report(left, right, args.tolerance, args.limit)

    if len(args.paths) != 1:
        p.error('give one directory to fingerprint, two directories to '
                'compare, or --compare with two .json files')

    fp = fingerprint(args.paths[0], args.frames)
    text = json.dumps(fp, indent=1, sort_keys=True)
    if args.output:
        with open(args.output, 'w') as h:
            h.write(text)
        print('%s: %d frames, %d arrays in the first, %d bytes'
              % (args.output, len(fp['frames']),
                 len(fp['frames'][0]['arrays']), len(text)))
    else:
        print(text)
    return 0


if __name__ == '__main__':
    sys.exit(main())
