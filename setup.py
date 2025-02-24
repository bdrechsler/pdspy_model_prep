from setuptools import setup, find_packages

setup(
    name="pdpsy_model_prep",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
        # e.g., 'requests>=2.25.1'
    ],
    author="Blake Drechsler",
    author_email="wbd814@gmail.com",
    description="Small package to help prep ALMA data to be fit with pdspy",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/my_project",  # Optional
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Choose appropriate license
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
