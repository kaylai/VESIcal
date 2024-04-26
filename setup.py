import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="VESIcal",
    version="1.2.10",
    author="Kayla Iacovino, Simon Matthews, Penny Wieser",
    author_email="kaylaiacovino@gmail.com",
    description=("A generalized python library for calculating and plotting various things "
                 "related to mixed volatile (H2O-CO2) solubility in silicate melts."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kaylai/VESIcal",
    packages=setuptools.find_packages(),
    install_requires=[
            'pandas',
            'numpy',
            'matplotlib',
            'scipy',
            'sympy',
            'openpyxl'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
