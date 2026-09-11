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
# Bumped whenever the comparison MATH changes, so two machines running
# different copies of this script say so instead of silently disagreeing.
#
# This is not hypothetical bookkeeping. The normalisation fix below changed the
# headline figure for one comparison from 3.608e-04 to 6.562e-08. A user ran the
# pre-fix copy against a post-fix analysis and got the old number with no
# indication why; the only tell was the absence of a line in the output. A tool
# whose whole purpose is comparing two machines should notice when the two
# machines are running different versions of it.
#
# The stamp is reported in the header of every comparison, so a pasted result
# identifies which maths produced it. It is deliberately NOT used to reject an
# older fingerprint file: v1 and v2 store the same five statistics, so an old
# .json is read correctly by a new script. Only the arithmetic moved, and that
# lives in the script, not the data - which is exactly why the output has to
# say so.
#
#   1  self-normalised statistics (the cancellation bug)
#   2  sorted_sum measured against abs_sorted_sum
COMPARISON_VERSION = 2

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

    out = {'source': os.path.abspath(directory),
           'comparison_version': COMPARISON_VERSION,
           'frames': []}
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


def _relative(a, b, scale=None):
    """Difference between two statistics, as a fraction of *scale*.

    WHY scale IS AN ARGUMENT, and why getting this wrong produced a false
    alarm on a 140-hour decision.

    `sorted_sum` of a near-zero-mean field is a catastrophic-cancellation
    quantity. Measured on the RMF antenna gate deck, frame 1, the cleaning
    potential phi: the sum of its 11989 values is 6.82e-10 while the sum of
    their magnitudes is 3.75e-06 - the sum is a 1.8e-04 residual of what went
    into it, 99.98% having cancelled.

    Normalising a difference in that sum BY THAT SUM therefore multiplies it by
    1/1.8e-04 = 5498. A real cross-build drift of 6.56e-08 was reported as
    3.608e-04 and tripped a 1e-10 tolerance, on a quantity whose absolute
    disagreement was 2.46e-13. Every one of the twelve worst entries was
    `sorted_sum`; not one `min`, `max` or `abs_sorted_sum` appeared, which is
    precisely the signature of cancellation rather than of a different answer.

    So a sum is measured against the magnitudes it was built from, not against
    the residue left after they cancelled. `min`, `max` and `abs_sorted_sum`
    do not cancel and keep their own scale.

    This is NOT a way of making the tool quieter: the absolute difference and
    the amplification are both reported, the cancellation-free statistics are
    still compared at full strength against the same tolerance, and a genuine
    difference moves min/max/abs_sorted_sum too. Verified against a case known
    to differ - the same deck at one rank and at two, which known-issues 15
    says gives a different answer - where the difference shows up at 1.9 on the
    sum AND 2.2e-02 of the magnitude, i.e. in both columns at once.

    WHERE THESE NUMBERS COME FROM, since this repository's first convention is
    that a claim names its evidence:

      3.19.6 magnitudes (6.8217e-10, 3.7506e-06)  measured here, from a run on
                                                  this machine
      the 3.25.5 side (2.461e-13, 6.562e-08, and
      that all twelve worst rows were sorted_sum) measured on the user's
                                                  cluster and relayed as a
                                                  pasted table; the .json
                                                  itself was never on this
                                                  machine
      1-rank vs 2-rank (1.9 / 0.826 / max)        measured here, both sides

    The distinction matters because an earlier commit message reported a
    'before and after' for the 3.25.5 comparison that was in fact produced by
    adding the relayed absolute differences to the local 3.19.6 fingerprint -
    a reconstruction that happens to reproduce the relayed numbers exactly, but
    which is not the same thing as running the tool on the real pair.
    """
    if a is None or b is None:
        return None
    if scale is None:
        scale = max(abs(a), abs(b))
    scale = abs(scale)
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
            # The field's own magnitude, used as the yardstick for the
            # cancelling statistic. See _relative().
            magnitude = max(ls.get('abs_sorted_sum') or 0.0,
                            rs.get('abs_sorted_sum') or 0.0)
            diffs = {k: _relative(ls[k], rs[k],
                                  magnitude if k == 'sorted_sum' else None)
                     for k in ('min', 'max', 'sorted_sum', 'abs_sorted_sum')}
            # Kept so the report can say how much the old normalisation would
            # have amplified this, rather than silently dropping the number.
            diffs['_sorted_sum_selfnorm'] = _relative(ls['sorted_sum'],
                                                      rs['sorted_sum'])
            diffs['_absdiff'] = abs((ls['sorted_sum'] or 0.0)
                                    - (rs['sorted_sum'] or 0.0))
            got = max((d for k, d in diffs.items()
                       if d is not None and not k.startswith('_')), default=0.0)
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
    print('comparison version %d' % COMPARISON_VERSION)
    print()

    # Deliberately NOT warning when a fingerprint carries an older stamp. The
    # file format has not changed - v1 and v2 store the same five statistics -
    # so an old .json is read correctly by a new script and needs no rewriting.
    # Only the COMPARISON math changed, and that lives in whichever copy of this
    # script is running, which is what the version line above reports.

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

    # Say which statistic produced that number and how big the disagreement
    # actually is. Without this the reader cannot tell a cancelling sum from a
    # moved extreme, and those mean opposite things.
    worst_diffs = fields[0][3]
    driving = max(((k, v) for k, v in worst_diffs.items()
                   if v is not None and not k.startswith('_')),
                  key=lambda kv: kv[1])[0]
    print('  driven by:      %s' % driving)
    if driving == 'sorted_sum':
        raw = worst_diffs.get('_sorted_sum_selfnorm')
        absd = worst_diffs.get('_absdiff')
        if absd is not None:
            print('  absolute:       %.3e' % absd)
        if raw and top > 0:
            print('  note:           this is a SUM over a field that largely')
            print('                  cancels. Measured against the sum itself')
            print('                  rather than the magnitudes behind it, the')
            print('                  same difference reads %.3e - %.0fx larger.'
                  % (raw, raw / top))
            print('                  The magnitudes are the honest yardstick.')
    moved = [k for k in ('min', 'max', 'abs_sorted_sum')
             if (worst_diffs.get(k) or 0.0) > tolerance]
    if driving == 'sorted_sum' and not moved:
        print('  and:            min, max and abs_sorted_sum all agree within')
        print('                  the tolerance. A different ANSWER moves those')
        print('                  too; only round-off moves the cancelling sum')
        print('                  alone.')
    if blocked:
        print('NOT COMPARABLE - see the refusals above. The number just')
        print('printed covers only the arrays that lined up, so it is not')
        print('evidence of anything.')
        return 2
    if top <= tolerance:
        print('within the %.1e tolerance: the two runs agree to round-off.'
              % tolerance)
        return 0
    print('ABOVE the %.1e tolerance.' % tolerance)
    if driving == 'sorted_sum' and not moved:
        # Do not decide this for the reader: give them the two regimes and the
        # number, and let them say which comparison they are making.
        print()
        print('WHICH REGIME IS THIS? The default tolerance assumes the two runs')
        print('should be bit-identical, which holds only for the same binary on')
        print('the same machine. Across a different compiler, BLAS or PETSc,')
        print('every operation reorders and round-off accumulates: on this')
        print('solver, ~1e-8 relative over a few hundred steps is ordinary.')
        print('So read it this way:')
        print('  same build, same machine  -> anything above 0 is a real bug')
        print('  different build or BLAS   -> ~1e-8 is expected; 1e-3 is not')
        print('Re-run with --tolerance if you are making the second comparison.')
        print('The absolute difference above is the number to judge, not this')
        print('ratio.')
        return 1
    print('These are not the same answer.')
    print('This is the signal the comparison exists to produce - do not')
    print('explain it away without finding the cause.')
    return 1


class _Refused(Exception):
    """Cannot compare at all - distinct from comparable-and-disagreeing.

    These must reach the caller as exit 2, not 1. A shell gate reads only the
    status, and "you typed the wrong path" and "the two builds disagree" call
    for opposite reactions; collapsing them into 1 makes a typo look like a
    finding and a finding look like a typo.
    """


def _load(path):
    """Read a fingerprint, refusing legibly rather than raising a traceback.

    The three ways this goes wrong in practice are a file that is not there
    (the reference was written on another machine and never copied over), a
    directory passed where a .json was expected, and a file that is not a
    fingerprint at all. A stack trace names none of them.
    """
    if os.path.isdir(path):
        raise _Refused(
            '%s is a directory. Either drop --compare to fingerprint two run\n'
            'directories directly, or pass the .json files written with -o.'
            % path)
    if not os.path.exists(path):
        raise _Refused(
            'no such fingerprint: %s\n'
            'Fingerprints are not produced by this flag - write one first on\n'
            'the machine that holds the run:\n'
            '    python3 scripts/run_fingerprint.py <results-dir> -o %s\n'
            'and copy it here if the two runs are on different machines.'
            % (path, os.path.basename(path)))
    try:
        with open(path) as handle:
            loaded = json.load(handle)
    except ValueError as exc:
        raise _Refused('%s is not valid JSON: %s' % (path, exc))
    if not isinstance(loaded, dict) or 'frames' not in loaded:
        raise _Refused(
            '%s is not a fingerprint - it has no "frames". Fingerprints are\n'
            'written by this script with -o; this looks like something else.'
            % path)
    return loaded


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
        return _report(_load(args.paths[0]), _load(args.paths[1]),
                       args.tolerance, args.limit)

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
    try:
        sys.exit(main())
    except _Refused as refusal:
        print(refusal, file=sys.stderr)
        sys.exit(2)
