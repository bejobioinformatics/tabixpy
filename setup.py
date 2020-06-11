import setuptools

with open("README.md", "rt") as fh:
    long_description = fh.read()

with open("VERSION", "rt") as fh:
    version = fh.read()


setuptools.setup(
    name="tabixpy", # Replace with your own username
    version=version,
    author="Saulo Aflitos",
    author_email="saulobejo@users.noreply.github.com",
    description="Tabix reader written 100% in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bejobioinformatics/tabixpy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
