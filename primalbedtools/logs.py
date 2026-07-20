"""Logging configuration for primalbedtools.

Library modules obtain a logger via ``get_logger(__name__)`` and never attach
handlers of their own. The CLI entrypoint calls ``configure_cli_logging`` to
route records to stderr, keeping stdout free for machine-readable output.

The package logger carries a NullHandler so library consumers who have not
configured logging see no output at all.
"""

import logging
import sys
from typing import Optional

LOGGER_NAME = "primalbedtools"

_package_logger = logging.getLogger(LOGGER_NAME)
_package_logger.addHandler(logging.NullHandler())

_cli_handler: Optional[logging.Handler] = None


def get_logger(name: str) -> logging.Logger:
    """Return a logger within the primalbedtools namespace.

    Args:
        name: Usually ``__name__``. Namespaced under the package logger if it
            is not already.

    Returns:
        logging.Logger: A child of the primalbedtools package logger.
    """
    if name == LOGGER_NAME or name.startswith(LOGGER_NAME + "."):
        return logging.getLogger(name)
    return logging.getLogger(f"{LOGGER_NAME}.{name}")


def configure_cli_logging(verbose: bool = False, quiet: bool = False) -> None:
    """Attach a stderr handler to the package logger. Safe to call repeatedly.

    Args:
        verbose: Emit debug-level detail.
        quiet: Suppress warnings, leaving only errors.
    """
    global _cli_handler

    if quiet:
        level = logging.ERROR
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.WARNING

    # The handler binds sys.stderr at construction, so it is built here rather
    # than at import time to respect any redirection in place when called.
    if _cli_handler is None:
        _cli_handler = logging.StreamHandler(sys.stderr)
        _cli_handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        _package_logger.addHandler(_cli_handler)

    _cli_handler.setLevel(level)
    _package_logger.setLevel(level)
    _package_logger.propagate = False  # the CLI owns its own output


def reset_cli_logging() -> None:
    """Undo configure_cli_logging, restoring library defaults."""
    global _cli_handler

    if _cli_handler is not None:
        _package_logger.removeHandler(_cli_handler)
        _cli_handler = None

    _package_logger.setLevel(logging.NOTSET)
    _package_logger.propagate = True
