# Spherical elastica: computer-assisted comparison

This directory contains the integer interval certificate used in Appendix A
of the research manuscript *A sharp elastic-energy inequality for spherical
Jordan curves*. It settles the remaining compact parameter range of the
exterior free-drop comparison.

The certificate supports one step of the geometric proof. It does not by
itself prove the contact-structure theorems or the whole spherical inequality.
The manuscript has undergone internal review, not independent peer review or
formal verification. No bibliographical priority claim is made.

## Quick replay

Requirements: Python 3.10 or later. Only the Python standard library is used.
Download or clone the repository, open a terminal in this directory, and run:

```sh
python audit_outside_certificate.py
```

Expected final fields:

```text
status: passed
coverage: complete prefix-free dyadic cover
leaf_count: 9367
max_depth: 18
arithmetic_checks: 12000
```

The recorded complete replay took about 232 seconds on the preparation
machine; runtime varies. The audit writes a fresh
`outside_one_cycle_certificate_audit.json` with the measured runtime.

## Files

| File | Purpose |
|---|---|
| `certify_outside_one_cycle.py` | Outward-rounded integer interval primitives and first certificate pass. |
| `refine_outside_certificate.py` | Refines the eight boxes left undecided by the first pass. |
| `audit_outside_certificate.py` | Checks the dyadic cover, replays all leaf inequalities and tests interval arithmetic against exact rational formulas. |
| `outside_one_cycle_certificate.json` | Original 48-cell pass; eight undecided boxes remain, so this file alone is not a complete certificate. |
| `outside_one_cycle_certificate_refined.json` | Complete certificate with 9,367 leaves and no undecided boxes. |
| `outside_one_cycle_certificate_audit.json` | Recorded successful replay of the complete certificate. |
| `SHA256SUMS` | SHA-256 hashes of the six proof/program/data files as uploaded. |

The audit deliberately reuses the generator's interval primitives. It is a
complete replay plus arithmetic checks, not a second independent implementation.

## Regenerate from scratch

Run these commands in order:

```sh
python certify_outside_one_cycle.py
python refine_outside_certificate.py
python audit_outside_certificate.py
```

The first command alone is insufficient. The refiner increases the integration
partition from 48 to 96 cells only for the eight undecided rectangles. Generated
JSON metadata contains runtime information, so regenerated file hashes can
differ even when the mathematical certificate is unchanged.

## Arithmetic and coverage

The fixed-point scale is `2**44`. Addition, subtraction, multiplication and
division are rounded outward; square-root bounds use integer square roots.
Pi is enclosed using rational alternating-series bounds in Machin's formula.
Floating point only proposes a trial cap parameter; its final rational value
is checked by integer interval arithmetic and need not be the optimal root.

Every leaf records its dyadic path, exact integer endpoints, rule and strict
margin. The audit checks prefix-freeness, the Kraft sum equal to one and endpoint
reconstruction, then recomputes every leaf decision. Counts are 96 small-pressure,
3 large-pressure, 93 elementary-energy, 1,037 integral-energy and 8,138
longitude leaves. The degenerate boundary's negative-phase angle pi is retained.

## 中文说明

本目录对应球面弹性能量研究稿附录 A 的计算机辅助比较。直接运行
`python audit_outside_certificate.py` 可复核完整的 9,367 个叶矩形及
12,000 项精确算术检查，无需第三方 Python 库。

若从头生成，必须按上面的三个命令依次运行；初次证书的 8 个未决矩形
须经过细化，不能把第一份 JSON 单独当成完整证明。复算器使用原来的
整数区间基本实现，不属于独立的第二套实现。这里的代码只承担论文中
指定的自由弧参数比较，几何归约仍依赖论文正文。
