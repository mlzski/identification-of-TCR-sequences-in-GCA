"""Install package."""
import io
import re

from pkg_resources import parse_requirements
from setuptools import find_packages, setup


def read_version(filepath: str) -> str:
    """Read the __version__ variable from the file.

    Args:
        filepath: probably the path to the root __init__.py

    Returns:
        the version
    """
    match = re.search(
        r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        io.open(filepath, encoding="utf_8_sig").read(),
    )
    if match is None:
        raise SystemExit("Version number not found.")
    return match.group(1)


# ease installation during development
vcs = re.compile(r"(git|svn|hg|bzr)\+")
try:
    with open("requirements.txt") as fp:
        VCS_REQUIREMENTS = [
            str(requirement)
            for requirement in parse_requirements(fp)
            if vcs.search(str(requirement))
        ]
except FileNotFoundError:
    # requires verbose flags to show
    print("requirements.txt not found.")
    VCS_REQUIREMENTS = []

# TODO: Update these values according to the name of the module.
setup(
    name="detectTCRs",
    version=read_version("detectTCRs/__init__.py"),  # single place for version
    description="Package for identification and analysis of TCR sequences associated with a disease.",
    long_description=open("README.md").read(),
    url="https://github.com/mlzski/identification-of-TCR-sequences-in-GCA",
    author="Anna Weber",
    author_email="a92.weber@t-online.de",
    maintainer="Michal Zulcinski",
    maintainer_email="m.zulcinski@leeds.ac.uk",
    # the following exclusion is to prevent shipping of tests.
    # if you do include them, add pytest to the required packages.
    packages=find_packages(".", exclude=["*tests*"]),
    install_requires=[
        # versions should be very loose here, just exclude unsuitable versions
        # because your dependencies also have dependencies and so on ...
        # being too strict here will make dependency resolution harder
    ],
)
