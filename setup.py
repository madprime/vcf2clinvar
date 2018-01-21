from setuptools import setup
from os import path

here = path.abspath(path.dirname(__file__))

setup(
    name='vcf2clinvar',
    version='0.1.5.0',
    description='Match a personal genome VCF datafile to ClinVar VCF',
    url='https://github.com/madprime/vcf2clinvar',
    author='Mad Price Ball',
    author_email='pgp-informatics@arvados.org',
    license='MIT',

    classifiers=[
        'Development Status :: 3 - Alpha',
    ],

    packages=['vcf2clinvar'],

    # Core dependencies should be listed here (will be installed by pip).
    install_requires=[],

    scripts=['clinvar_report.py'],

)
