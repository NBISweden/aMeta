# Contributor guide

Thank you for your interest in improving this project. This project is
open-source under the [MIT license] and welcomes contributions in the
form of bug reports, feature requests, and pull requests.

Here is a list of important resources for contributors:

- [Source Code]
- [Issue Tracker]

[mit license]: https://opensource.org/licenses/MIT
[source code]: https://github.com/NBISweden/aMeta
[issue tracker]: https://github.com/NBISweden/aMeta/issues

## Development environment

Follow the installation instructions in the README. Checkout a
development branch from main and apply changes. Use the nox tests to
test the integrity of the workflow.

## Nox tests

In addition to the regular tests, there is also the option to run
tests based on the [nox](https://nox.thea.codes) framework. These
tests are mainly aimed at developers.

### Requirements

To run nox tests you need to install a number of requirements. See the
[nox homepage](https://nox.thea.codes) for details; a short summary
follows here.

Install nox, either with [pipx](https://pipx.pypa.io/stable/)

    pipx install nox

or [pip](https://pip.pypa.io/en/stable/)

    python3 -m pip install nox

### Nox sessions

The test file `.test/noxfile.py` defines parametrized test sessions
that you can list with

    nox --list

nox will test different combinations of Python and Snakemake and run
[pytest](https://docs.pytest.org/en/8.2.x/) tests in isolated test
directories. You can [specify parametrized
sessions](https://nox.thea.codes/en/stable/usage.html#specifying-parametrized-sessions)
by passing the session name:

    nox --session "snakemake(python='3.11', snakemake='7.32.4')"

To reuse nox sessions and avoid reinstallation of dependencies, use
the `-R` option.
