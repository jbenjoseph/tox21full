from setuptools import setup, find_packages

with open("requirements.txt") as fd:
    install_requires = fd.read().splitlines()

setup(
    name="tox21full",
    version="0.1.0",
    description="Creates the complete Tox21 dataset",
    long_description=open("README.rst").read(),
    keywords="machine_learning artificial_intelligence",
    author="JJ Ben-Joseph",
    author_email="jj@memoriesofzion.com",
    python_requires=">=3.8",
    url="https://github.com/etz4ai/tox21full",
    license="Apache",
    classifiers=[
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3 :: Only",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    install_requires=install_requires,
    entry_points={"console_scripts": ["tox21full = tox21full.__main__:main"]},
)
