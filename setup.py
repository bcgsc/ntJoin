import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="minimizer_assemble",
    version="0.1dev",
    author="Lauren Coombe",
    author_email="lcoombe@bcgsc.ca",
    description="Assembly using minimizers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bcgsc/minimizer-assembly-scaffolder",
    license="GPLv3",
    python_requires=">=3",
    install_requires=["networkx", "pybedtools"],
    packages = ["minimizer-assembly-scaffolder"],
    scripts = ["minimizerize.py", "minimizer_assemble.py", "hash.py", "read_fasta.py"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)