import os
import re
from ftplib import FTP

DIR_CLINVAR_VCF_B37 = 'pub/clinvar/vcf_GRCh37'
DIR_CLINVAR_VCF_B38 = 'pub/clinvar/vcf_GRCh38'


def nav_to_vcf_dir(ftp, build='b37'):
    # Navigate to appropriate directory.
    if build == 'b37':
        ftp.cwd(DIR_CLINVAR_VCF_B37)
    elif build == 'b38':
        ftp.cwd(DIR_CLINVAR_VCF_B38)
    else:
        raise IOError("Genome build not recognized.")


def latest_vcf_filename(build='b37'):
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    nav_to_vcf_dir(ftp, build=build)
    clinvar_datestamped = [f for f in ftp.nlst() if
                           re.match('^clinvar_[0-9]{8}.vcf.gz$', f)]
    if len(clinvar_datestamped) == 1:
        return clinvar_datestamped[0]
    raise IOError("Unable to determine the most recent ClinVar VCF file on " +
                  "NCBI's FTP site.")


def get_latest_vcf_file(target_dir, build='b37', overwrite=False):
    latest_filename = latest_vcf_filename()
    target_filename = os.path.join(target_dir, latest_filename)

    # Abort if we believe we already have this file and overwrite is False.
    if not overwrite and os.path.exists(target_filename):
        return

    # Get the file.
    with open(target_filename, 'w') as fh:
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        nav_to_vcf_dir(ftp, build=build)
        ftp.retrbinary('RETR %s' % latest_filename, fh.write)
