"""Performance monitoring for SSWD-EvoEpi simulations.

Provides component-level timing instrumentation that can be enabled
with a single flag. Designed to be non-invasive — zero overhead when
disabled.

Usage:
    from sswd_evoepi.perf import PerfMonitor

    perf = PerfMonitor(enabled=True)
    
    with perf.track("disease"):
        daily_disease_update(...)
    
    with perf.track("spawning"):
        spawning_step(...)
    
    perf.report()  # Print component breakdown
"""

import time
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass
class ComponentStats:
    """Timing statistics for a single component."""
    total_time: float = 0.0
    call_count: int = 0
    min_time: float = float('inf')
    max_time: float = 0.0

    @property
    def mean_time(self) -> float:
        return self.total_time / self.call_count if self.call_count > 0 else 0.0

    @property
    def pct_of(self) -> float:
        """Set externally by report()."""
        return 0.0


class PerfMonitor:
    """Lightweight component-level performance monitor.
    
    When disabled, all methods are no-ops with zero overhead.
    When enabled, tracks wall-clock time per named component.
    """

    def __init__(self, enabled: bool = False):
        self.enabled = enabled
        self._stats: Dict[str, ComponentStats] = defaultdict(ComponentStats)
        self._start_time: Optional[float] = None
        self._total_time: float = 0.0

    def start(self) -> None:
        """Mark the start of the monitored region."""
        if self.enabled:
            self._start_time = time.perf_counter()

    def stop(self) -> None:
        """Mark the end of the monitored region."""
        if self.enabled and self._start_time is not None:
            self._total_time = time.perf_counter() - self._start_time

    @contextmanager
    def track(self, component: str):
        """Context manager to time a named component.
        
        Usage:
            with perf.track("disease"):
                do_disease_stuff()
        """
        if not self.enabled:
            yield
            return

        t0 = time.perf_counter()
        yield
        elapsed = time.perf_counter() - t0

        stats = self._stats[component]
        stats.total_time += elapsed
        stats.call_count += 1
        stats.min_time = min(stats.min_time, elapsed)
        stats.max_time = max(stats.max_time, elapsed)

    def record(self, component: str, elapsed: float) -> None:
        """Manually record a timing measurement."""
        if not self.enabled:
            return
        stats = self._stats[component]
        stats.total_time += elapsed
        stats.call_count += 1
        stats.min_time = min(stats.min_time, elapsed)
        stats.max_time = max(stats.max_time, elapsed)

    def get_stats(self) -> Dict[str, ComponentStats]:
        """Return raw component statistics."""
        return dict(self._stats)

    def summary(self) -> dict:
        """Return a summary dict suitable for JSON serialization."""
        total = self._total_time or sum(s.total_time for s in self._stats.values())
        result = {}
        for name, stats in sorted(self._stats.items(), key=lambda x: -x[1].total_time):
            pct = (stats.total_time / total * 100) if total > 0 else 0
            result[name] = {
                'total_s': round(stats.total_time, 4),
                'calls': stats.call_count,
                'mean_ms': round(stats.mean_time * 1000, 3),
                'pct': round(pct, 1),
            }
        result['_total_s'] = round(total, 4)
        return result

    def report(self, title: str = "Performance Breakdown") -> str:
        """Generate a human-readable performance report."""
        total = self._total_time or sum(s.total_time for s in self._stats.values())
        
        lines = [
            f"\n{'='*60}",
            f" {title}",
            f"{'='*60}",
            f"{'Component':<25} {'Total (s)':>10} {'Calls':>8} {'Mean (ms)':>10} {'%':>6}",
            f"{'-'*25} {'-'*10} {'-'*8} {'-'*10} {'-'*6}",
        ]

        for name, stats in sorted(self._stats.items(), key=lambda x: -x[1].total_time):
            pct = (stats.total_time / total * 100) if total > 0 else 0
            lines.append(
                f"{name:<25} {stats.total_time:>10.4f} {stats.call_count:>8} "
                f"{stats.mean_time*1000:>10.3f} {pct:>5.1f}%"
            )

        lines.append(f"{'-'*25} {'-'*10} {'-'*8} {'-'*10} {'-'*6}")
        lines.append(f"{'TOTAL':<25} {total:>10.4f}")
        lines.append(f"{'='*60}\n")

        return '\n'.join(lines)

    def reset(self) -> None:
        """Clear all accumulated statistics."""
        self._stats.clear()
        self._start_time = None
        self._total_time = 0.0
