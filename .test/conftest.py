import pytest

def pytest_addoption(parser):
    parser.addoption("--num-cores", action="store", default=4,
                     type=int, help="Number of cores to use")

@pytest.fixture
def num_cores(request):
    ncores = request.config.getoption("--num-cores")
    if ncores < 1:
        raise ValueError("Number of cores must be greater than 0")
    return ncores
