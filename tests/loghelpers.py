import contextlib
import logging


class _ListHandler(logging.Handler):
    """Collects records rather than emitting them."""

    def __init__(self):
        super().__init__()
        self.records = []

    def emit(self, record):
        self.records.append(record)


@contextlib.contextmanager
def capture_logs(name="primalbedtools", level=logging.DEBUG):
    """Capture records from a primalbedtools logger.

    unittest gained assertNoLogs in 3.10, but this project supports 3.9, so this
    is the way to assert that nothing was logged.
    """
    logger = logging.getLogger(name)
    handler = _ListHandler()

    old_level = logger.level
    old_propagate = logger.propagate
    old_disabled = logger.disabled

    logger.setLevel(level)
    logger.propagate = False
    logger.disabled = False
    logger.addHandler(handler)
    try:
        yield handler.records
    finally:
        logger.removeHandler(handler)
        logger.setLevel(old_level)
        logger.propagate = old_propagate
        logger.disabled = old_disabled
