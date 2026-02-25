from setuptools import setup

setup(
    name='sidle-standalone',
    version='1.0.0',
    description='Multi-region 16S reconstruction without QIIME2',
    author='Modified from q2-sidle by Justine Debelius',
    license='BSD-3-Clause',
    packages=['sidle_standalone'],
    package_dir={'sidle_standalone': '.'},
    install_requires=[
        'numpy>=1.19.0',
        'pandas>=1.1.0',
        'biopython>=1.78',
        'scikit-bio>=0.5.6',
    ],
    entry_points={
        'console_scripts': [
            'sidle=sidle_standalone.cli:main',
        ],
    },
    python_requires='>=3.7',
)
