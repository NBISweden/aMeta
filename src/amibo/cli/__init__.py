import os

ROOT = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
TESTDIR = os.path.join(ROOT, ".test")
SNAKEFILE = os.path.join(ROOT, "workflow", "Snakefile")
