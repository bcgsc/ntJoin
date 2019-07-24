import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ntJoin",
    version="0.1dev",
    author="Lauren Coombe",
    author_email="lcoombe@bcgsc.ca",
    description="Assembly using minimizers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bcgsc/ntJoin",
    license="GPLv3",
    python_requires=">=3",
    install_requires=["python-igraph", "pybedtools"],
    scripts = ["bin/ntjoin_assemble.py", "bin/read_fasta.py"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)