"""Run logging: timestamped progress to a file in the repository, with an ETA.

Long runs used to write to a session-scoped scratch directory and go silent for tens of
minutes at a time, which makes it impossible to tell a slow job from a dead one. A job now
writes to ``logs/<name>.log`` inside the repository and emits a line for every completed
point, carrying elapsed time and a projected finish.

    from berrycasscf.runlog import JobLog

    log = JobLog("butadiene_scan", total=65)
    scan_gap(..., progress=log)
    log.done()

``JobLog`` is itself the progress callback, so it drops into the ``progress=`` argument the
drivers already take.
"""

from __future__ import annotations

import os
import sys
import time
from datetime import datetime, timedelta

DEFAULT_LOG_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "logs"
)


def _hms(seconds: float) -> str:
    return str(timedelta(seconds=int(max(seconds, 0))))


class JobLog:
    """A progress callback that timestamps to stdout and to ``logs/<name>.log``.

    ``total``, when known, turns each completed step into an ETA. Nothing here fails a run:
    if the log file cannot be opened the job still proceeds, writing to stdout only.
    """

    def __init__(self, name: str, total: int | None = None, log_dir: str | None = None,
                 echo: bool = True):
        self.name = name
        self.total = total
        self.echo = echo
        self.started = time.time()
        self.count = 0
        self._fh = None
        directory = log_dir or DEFAULT_LOG_DIR
        try:
            os.makedirs(directory, exist_ok=True)
            self.path = os.path.join(directory, f"{name}.log")
            self._fh = open(self.path, "a", buffering=1)      # line buffered
        except OSError as exc:                                 # pragma: no cover - rare
            self.path = None
            print(f"[{name}] could not open log file ({exc}); logging to stdout only")
        self._emit(f"=== {name} started {datetime.now():%Y-%m-%d %H:%M:%S}"
                   + (f", {total} steps expected" if total else "") + " ===")

    # -- the progress callback -----------------------------------------------------
    def __call__(self, message: str) -> None:
        """Log one line. Counts as a completed step when ``total`` is known."""
        if self.total:
            self.count += 1
            self._emit(f"[{self.count}/{self.total}] {message}{self._eta()}")
        else:
            self._emit(str(message))

    def note(self, message: str) -> None:
        """Log without advancing the step counter."""
        self._emit(str(message))

    def done(self, message: str = "finished") -> None:
        self._emit(f"=== {self.name} {message} after {_hms(time.time()-self.started)} ===")
        if self._fh:
            self._fh.close()
            self._fh = None

    # -- internals -----------------------------------------------------------------
    def _eta(self) -> str:
        if not self.total or self.count == 0:
            return ""
        elapsed = time.time() - self.started
        per = elapsed / self.count
        remaining = per * (self.total - self.count)
        if self.count >= self.total:
            return f"  [{_hms(elapsed)} elapsed]"
        finish = datetime.now() + timedelta(seconds=remaining)
        return (f"  [{_hms(elapsed)} elapsed, {per:.1f} s/step, "
                f"~{_hms(remaining)} left, ETA {finish:%H:%M}]")

    def _emit(self, line: str) -> None:
        stamped = f"{datetime.now():%H:%M:%S} {line}"
        if self._fh:
            self._fh.write(stamped + "\n")
        if self.echo:
            print(stamped, flush=True)

    def __enter__(self) -> "JobLog":
        return self

    def __exit__(self, *exc) -> None:
        self.done("failed" if exc[0] else "finished")
