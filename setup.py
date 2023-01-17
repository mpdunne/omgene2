from setuptools import setup, find_packages

with open("requirements.txt", "r") as file_handler:
    requirements = file_handler.readlines()

setup(
    name="omgene",
    version="0.0.0",
    include_package_data=True,
    install_requires=requirements,
    packages=find_packages(include=["omgene", "omgene.*"]),
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
