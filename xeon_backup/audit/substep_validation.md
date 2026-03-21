# Substep Validation: substeps=6 vs substeps=24

**Date:** 2026-02-26 09:25

## Methodology

Single-seed comparisons are meaningless because changing substeps
changes the number of RNG draws per day, shifting the entire random
sequence. Instead, we compare **distributions** across many seeds.

- **Seeds:** 0–9 (10 per setting)
- **Network:** 5 nodes, K scaled to ~2000 total
- **Duration:** 5 years, disease introduced at year 1
- **Total runs:** 20

## Results Summary

| Metric | substeps=6 (mean±std) | substeps=24 (mean±std) | Welch p | M-W p | Cohen's d | Significant? |
|--------|-----------------------|------------------------|---------|-------|-----------|-------------|
| final_pop | 1637.6 ± 68.1 | 1374.4 ± 112.8 | 0.0000 | 0.0002 | 2.825 | Yes |
| crash_pct | 18.5 ± 3.0 | 31.2 ± 5.6 | 0.0000 | 0.0002 | -2.800 | Yes |
| disease_deaths | 1603.1 ± 42.2 | 1924.0 ± 82.9 | 0.0000 | 0.0002 | -4.878 | Yes |
| min_pop | 1627.9 ± 60.6 | 1374.4 ± 112.8 | 0.0000 | 0.0002 | 2.800 | Yes |

## Runtime

- **substeps=6:** 8.2s mean
- **substeps=24:** 12.3s mean
- **Speedup:** 1.51×

## Raw Data

### substeps=6

| Seed | Final Pop | Crash % | Disease Deaths | Min Pop | Runtime (s) |
|------|-----------|---------|----------------|---------|-------------|
| 0 | 1726 | 18.2% | 1647 | 1635 | 8.4 |
| 1 | 1584 | 20.7% | 1632 | 1584 | 8.3 |
| 2 | 1638 | 18.3% | 1623 | 1632 | 8.4 |
| 3 | 1550 | 22.4% | 1573 | 1550 | 8.0 |
| 4 | 1742 | 12.8% | 1658 | 1742 | 8.0 |
| 5 | 1605 | 19.7% | 1616 | 1605 | 8.1 |
| 6 | 1615 | 19.2% | 1568 | 1615 | 8.1 |
| 7 | 1569 | 21.5% | 1606 | 1569 | 8.2 |
| 8 | 1716 | 14.1% | 1591 | 1716 | 8.1 |
| 9 | 1631 | 18.4% | 1517 | 1631 | 8.2 |

### substeps=24

| Seed | Final Pop | Crash % | Disease Deaths | Min Pop | Runtime (s) |
|------|-----------|---------|----------------|---------|-------------|
| 0 | 1414 | 29.2% | 2041 | 1414 | 12.3 |
| 1 | 1394 | 30.2% | 1800 | 1394 | 11.9 |
| 2 | 1414 | 29.2% | 1836 | 1414 | 12.4 |
| 3 | 1394 | 30.2% | 1804 | 1394 | 12.3 |
| 4 | 1472 | 26.3% | 1977 | 1472 | 12.3 |
| 5 | 1371 | 31.4% | 1970 | 1371 | 12.4 |
| 6 | 1546 | 22.6% | 1987 | 1546 | 12.6 |
| 7 | 1307 | 34.6% | 1925 | 1307 | 12.6 |
| 8 | 1124 | 43.7% | 1964 | 1124 | 12.3 |
| 9 | 1308 | 34.5% | 1936 | 1308 | 12.4 |

## Per-Seed Comparison (node-level final populations)

**Seed 0:**
- substeps=6:  [541, 191, 258, 341, 395] (total=1726)
- substeps=24: [540, 180, 84, 214, 396] (total=1414)

**Seed 1:**
- substeps=6:  [534, 216, 267, 209, 358] (total=1584)
- substeps=24: [489, 182, 103, 223, 397] (total=1394)

**Seed 2:**
- substeps=6:  [539, 205, 312, 252, 330] (total=1638)
- substeps=24: [512, 191, 194, 126, 391] (total=1414)

**Seed 3:**
- substeps=6:  [505, 202, 266, 184, 393] (total=1550)
- substeps=24: [514, 181, 128, 174, 397] (total=1394)

**Seed 4:**
- substeps=6:  [522, 198, 285, 339, 398] (total=1742)
- substeps=24: [513, 191, 273, 154, 341] (total=1472)

**Seed 5:**
- substeps=6:  [509, 168, 325, 206, 397] (total=1605)
- substeps=24: [502, 172, 182, 120, 395] (total=1371)

**Seed 6:**
- substeps=6:  [536, 197, 276, 210, 396] (total=1615)
- substeps=24: [505, 205, 239, 204, 393] (total=1546)

**Seed 7:**
- substeps=6:  [516, 198, 291, 167, 397] (total=1569)
- substeps=24: [505, 187, 257, 199, 159] (total=1307)

**Seed 8:**
- substeps=6:  [534, 189, 257, 342, 394] (total=1716)
- substeps=24: [498, 162, 184, 172, 108] (total=1124)

**Seed 9:**
- substeps=6:  [534, 194, 256, 254, 393] (total=1631)
- substeps=24: [535, 206, 217, 242, 108] (total=1308)

## Verdict

**NO — substeps=6 produces significantly different dynamics**

Statistically significant differences with large effect sizes were found.
substeps=24 should be retained for accurate dynamics.
