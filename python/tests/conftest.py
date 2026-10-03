"""Shared pytest configuration.

Tests that take more than a few seconds are marked ``@pytest.mark.slow`` and
skipped by default.  Run everything with ``pytest --runslow`` (about 15 minutes),
or only the slow ones with ``pytest --runslow -m slow``.
"""
import pytest


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False,
                     help="also run tests marked slow (full pipeline / high resolution)")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: takes more than a few seconds; run with --runslow")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        return
    skip = pytest.mark.skip(reason="slow: use --runslow")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)
