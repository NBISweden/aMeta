#!/usr/bin/env python3
import os

import pytest


def pytest_configure(config):
    pytest.dname = os.path.dirname(__file__)
    pytest.project = os.path.dirname(pytest.dname)


@pytest.fixture(autouse=False)
def cd_tmp_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
