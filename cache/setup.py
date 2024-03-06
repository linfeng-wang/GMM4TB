from setuptools import find_packages, setup

with open("app/README.md", "r") as f:
    long_description = f.read()

setup(
    name="GMM4TB",
    version="0.0.10",
    description="A model for strain calling and drug resistance profiling",
    package_dir={"": "."},
    packages=find_packages(where="."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/linfeng-wang/GMM4TB",
    author="linfeng-wang",
    author_email="linfeng.wang@lshtm.ac.uk",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    install_requires=["bson >= 0.5.10"],
    extras_require={
        "dev": ["pytest>=7.0", "twine>=4.0.2"],
    },
    python_requires=">=3.10",
)