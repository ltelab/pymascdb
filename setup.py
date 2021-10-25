import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mascdb",
    version="0.1.0",
    author="Grazioli",
    author_email="jacopo.grazioli@epfl.ch",
    description="A package to read and manipulate data of MASCdb database (snowflake images and retrievals)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jacgraz/pymascdb",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
