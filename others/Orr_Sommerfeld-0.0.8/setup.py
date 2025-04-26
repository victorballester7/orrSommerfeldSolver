import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Orr_Sommerfeld", #lib name (debe coincidir con folder_name)
    version="0.0.8",
    author="Fernando Scarafia",
    author_email="fernando.scarafia@ib.edu.ar",
    description="La presente libreria permite analizar Orr Sommerfeld temporal y espacial",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab/publish",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License"
    ),
)
