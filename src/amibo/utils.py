import contextlib
import os


@contextlib.contextmanager
def cd(path):
    """
    Convenience function to change the current working directory in a context manager.
    """
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)
