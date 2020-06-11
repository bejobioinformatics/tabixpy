import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tabixpy", # Replace with your own username
    version="1",
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
