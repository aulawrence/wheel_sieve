import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="wheel_sieve",
    version="0.0.1",
    author="Lawrence Au",
    author_email="aulawrence@gmail.com",
    description="Wheel Sieve and More",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/aulawrence/wheel_sieve",
    packages=["wheel_sieve"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.5",
    install_requires=["numpy", "gmpy2"],
)
