from setuptools import setup, find_packages

setup(
    name="disto",
    version="0.1",
    packages=find_packages(),  # This will find both src and disto_cli packages
    package_dir={
        'src': 'src',
        'disto_cli': 'disto_cli'
    },
    entry_points={
        'console_scripts': [
            'disto=disto_cli.cli:main',
        ],
    },
)

