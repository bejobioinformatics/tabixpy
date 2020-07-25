import os
import sys
import setuptools

sys.path.insert(0, '.')

import tabixpy

setuptools.setup(
    name                          = tabixpy.__package_name__,
    version                       = tabixpy.__version__,
    author                        = tabixpy.__author__,
    author_email                  = tabixpy.__author_email__,
    description                   = tabixpy.__description__,
    long_description              = tabixpy.__long_description__,
    long_description_content_type = "text/markdown",
    # license                       = tabixpy.__license__,
    url                           = tabixpy.__url__,
    packages                      = setuptools.find_packages(),
    keywords                      = tabixpy.___keywords__,
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points={
        "console_scripts": [
            "tabixpy = tabixpy._main:main",
        ]
    },
    include_package_data=True,
)
