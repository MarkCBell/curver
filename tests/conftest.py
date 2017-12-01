
import pytest

def pytest_addoption(parser):
    parser.addoption('--runslow', action='store_true', help='run slow tests')

def pytest_collection_modifyitems(config, items):
    if not config.getoption('--runslow'):
        skip_slow = pytest.mark.skip(reason='need --runslow option to run')
        for item in items:
            if 'slow' in item.keywords:
                item.add_marker(skip_slow)

