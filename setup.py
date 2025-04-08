from setuptools import setup, find_packages
from pathlib import Path

# Read requirements from requirements.txt
def read_requirements():
    requirements_file = Path(__file__).parent / "requirements.txt"
    if requirements_file.exists():
        with open(requirements_file) as f:
            return f.read().splitlines()
    return []

setup(
    name="pyROX",
    version="1.0.0",
    author="Sam de Regt",
    author_email="regt@strw.leidenuniv.nl",
    description="pyROX: Rapid Opacity X-sections",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/samderegt/pyROX",  # Update with your repository URL
    packages=find_packages(),  # Automatically find package directories
    install_requires=read_requirements(),  # Read dependencies from requirements.txt
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    entry_points={
        "console_scripts": [
            "pyROX=pyROX.main:main", # Execute the main function in pyROX/main.py
        ],
    },
)