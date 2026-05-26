"""Drive Protenix v2 inference directly via the bundled Modal app.

This is a thin convenience wrapper for the user task "run Protenix
alone" — the same work happens inside ``run_complex.py
--models protenix``. The Modal app at
``scripts/modal/modal_protenix_avh.py`` is the actual workhorse.

Why a separate file: Protenix v2 (ByteDance, released 2026-04-08) is
the newest of the three predictors and the most likely to be swapped
in / out as the weight access situation evolves. Today the public PyPI
package ships ``protenix>=2.0.0`` but the ``protenix-v2`` checkpoint
itself is gated behind ByteDance's internal review (see
https://github.com/bytedance/Protenix/issues/295). The Modal app
defaults to ``protenix_base_default_v1.0.0`` (the AlphaFold-3-equivalent
public weight) and switches to v2 when ``PROTENIX_MODEL_NAME=protenix-v2``
is set.

Usage::

    mise run rerun:protenix
    # equivalently:
    uv run python scripts/data/run_protenix.py

    # Flip to v2 weights (once accessible)
    PROTENIX_MODEL_NAME=protenix-v2 mise run rerun:protenix
"""

from __future__ import annotations

import sys

from scripts.data.run_complex import main as complex_main


def main(argv: list[str] | None = None) -> None:
    """Run only the Protenix branch of the complex-rerun panel.

    Forwards every flag straight through to ``run_complex.main``
    (so ``--help``, ``--limit``, ``--download-only``, ``--launch-only``
    all work) and only injects ``--models protenix`` when the user
    hasn't passed ``--models`` themselves.
    """
    forwarded = list(argv) if argv is not None else list(sys.argv[1:])
    if "--models" not in forwarded:
        forwarded.extend(["--models", "protenix"])
    complex_main(forwarded)


if __name__ == "__main__":
    main()
