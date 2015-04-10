"""
Fuctions for automatically finding and retrieving recent ClinVar VCF file.

Notes:
- 'build' is explicitly required by these functions to avoid potential
    accidental build mismatches during processing.
- The most recent file is stored with a filename constructed by concatenating
    the genome build with the ClinVar file's original name. Example: if the
    ClinVar filename is 'clinvar_20150330.vcf.gz' and the build is 'b37', then
    the locally stored file is named 'b37_clinvar_20150330.vcf.gz'.
"""
import os
import re
from ftplib import FTP

DIR_CLINVAR_VCF_B37 = 'pub/clinvar/vcf_GRCh37'
DIR_CLINVAR_VCF_B38 = 'pub/clinvar/vcf_GRCh38'


def nav_to_vcf_dir(ftp, build):
    """
    Navigate an open ftplib.FTP to appropriate directory for ClinVar VCF files.

    Args:
        ftp:   (type: ftplib.FTP) an open connection to ftp.ncbi.nlm.nih.gov
        build: (type: string) genome build, either 'b37' or 'b38'
    """
    if build == 'b37':
        ftp.cwd(DIR_CLINVAR_VCF_B37)
    elif build == 'b38':
        ftp.cwd(DIR_CLINVAR_VCF_B38)
    else:
        raise IOError("Genome build not recognized.")


def latest_vcf_filename(build):
    """
    Determine the filename for the most recent comprehensive ClinVar VCF.

    Args:
        build: (type: string) genome build, either 'b37' or 'b38'

    Returns:
        (type: string) Filename of the most recent comprehensive ClinVar VCF.
    """
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    nav_to_vcf_dir(ftp, build=build)
    clinvar_datestamped = [f for f in ftp.nlst() if
                           re.match('^clinvar_[0-9]{8}.vcf.gz$', f)]
    if len(clinvar_datestamped) == 1:
        return clinvar_datestamped[0]
    raise IOError("Unable to determine the most recent ClinVar VCF file on " +
                  "NCBI's FTP site.")


def get_latest_vcf_file(target_dir, build, overwrite=False):
    """
    Download most recent ClinVar VCF to target local dir, if not already there.

    Args:
        target_dir: (type: string) Relative path to target local directory
        build:      (type: string) genome build, either 'b37' or 'b38'
        overwrite:  (type: boolean; optional, default=False) Whether to
                    a local ClinVarFile that appears to be current based on the
                    filename.

    Returns:
        (type: string) Local path to the up-to-date comprehensive ClinVar VCF.
    """
    latest_filename = latest_vcf_filename(build)
    target_filename = os.path.join(target_dir, build + '_' + latest_filename)

    # Abort if we believe we already have this file and overwrite is False.
    if not overwrite and os.path.exists(target_filename):
        return target_filename

    # Get the file.
    with open(target_filename, 'w') as fh:
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        nav_to_vcf_dir(ftp, build=build)
        ftp.retrbinary('RETR %s' % latest_filename, fh.write)

    return target_filename
