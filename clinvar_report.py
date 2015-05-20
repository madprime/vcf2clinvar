#!/usr/bin/python
#

from __future__ import unicode_literals
from argparse import ArgumentParser
import bz2
import csv
import gzip
import json
import os
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
import sys

from vcf2clinvar import match_to_clinvar
from vcf2clinvar.clinvar_update import get_latest_vcf_file

# Don't stack dump on keyboard ctrl-c or on premature
# termination of output stream (say from piping output
# through head).
#
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)


def json_report(input_genome_file, input_clinvar_file, build, notes, version):
    json_report = {}
    json_report["vcf2clinvar-version"] = version
    json_report["notes"] = notes
    json_report["genome-build"] = build
    json_report["variants"] = []
    matching = match_to_clinvar(input_genome_file, input_clinvar_file)
    for genome_vcf_line, clinvar_allele, zygosity in matching:
        for record in clinvar_allele.records:
            data = {
                'chrom': genome_vcf_line.chrom,
                'pos': genome_vcf_line.start,
                'ref-allele': genome_vcf_line.ref_allele,
                # Note that var_allele CAN be the same as the ref_allele, if
                # the reference variant is considered pathogenic! Factor V
                # Leiden is a famous example.
                'var-allele': clinvar_allele.sequence,
                'zygosity': zygosity,
            }
            json_report['variants'].append(data)
    print(json.dumps(json_report))


def csv_report(input_genome_file, input_clinvar_file, build, version):
    print('##genome-build=%s' % build)
    print('##vcf2clinvar-version=%s' % version)
    csv_out = csv.writer(sys.stdout, lineterminator=os.linesep)
    header = ('Chromosome', 'Position', 'Reference allele', 'Variant allele',
              'Name', 'Significance', 'Frequency', 'Zygosity',
              'ClinVar accession URL')
    csv_out.writerow(header)
    matching = match_to_clinvar(input_genome_file, input_clinvar_file)
    for genome_vcf_line, clinvar_allele, zygosity in matching:
        chrom = genome_vcf_line.chrom
        pos = genome_vcf_line.start
        ref_allele = genome_vcf_line.ref_allele
        alt_allele = clinvar_allele.sequence
        name = ':'.join([r.dbn for r in clinvar_allele.records])
        allele_freq = clinvar_allele.frequency
        zygosity = zygosity
        for record in clinvar_allele.records:
            name = record.dbn
            clnsig = record.sig
            url = url = 'http://www.ncbi.nlm.nih.gov/clinvar/' + record.acc
            data = (chrom, pos, ref_allele, alt_allele, name, clnsig,
                    allele_freq, zygosity, url)
            csv_out.writerow(data)


def main():
    """
    Parse command line argument and
    output appropriate file type (csv or JSON)
    """
    parser = ArgumentParser()

    parser.add_argument(
        "-c", "--clinvarfile", dest="clinvarfile",
        help="ClinVar VCF file (either this or -C must be specified)",
        metavar="CLINVARFILE")
    parser.add_argument(
        "-C", "--clinvardir", dest="clinvardir",
        help="ClinVar VCF directory (either this or -c must be specified). " +
        "This option will use vcf2clinvar.clinvar_update to automatically " +
        "check and import the most recent ClinVar file to this directory.",
        metavar="CLINVARDIR")
    parser.add_argument(
        "-i", "--input", dest="inputfile",
        help="Input VCF file ['.vcf', '.vcf.gz', '.vcf.bz2']. " +
        "Uncompressed genome data is also accepted via stdin.",
        metavar="INPUT")
    parser.add_argument(
        "-t", "--type", dest="type", default='csv',
        help="Output report type ('csv' or 'json'). Defaults to csv. " +
        "CSV Report: Reports all genome variants matching ClinVar records, " +
        "and some summary ClinVar data from these records. Header lines " +
        "with metadata begin with '##'.\n" +
        "JSON Report: Reports genome variants matching ClinVar records " +
        "(no record information is included).",
        metavar="TYPE")
    parser.add_argument(
        "-n", "--notes", dest="notes",
        help="Notes (JSON format) to include in report. (JSON report only)",
        metavar="NOTES")
    parser.add_argument(
        "-g", "--genome-build", dest="build",
        help="Genome build to include in report ('b37' or 'b38').",
        metavar="GENOMEBUILD")
    options = parser.parse_args()

    version = os.popen("python setup.py --version").read().strip()

    if not sys.stdin.isatty():
        input_genome_file = sys.stdin
    elif options.inputfile:
        if options.inputfile.endswith('.vcf'):
            input_genome_file = open(options.inputfile)
        elif options.inputfile.endswith('.vcf.gz'):
            input_genome_file = gzip.open(options.inputfile)
        elif options.inputfile.endswith('.vcf.bz2'):
            input_genome_file = bz2.BZ2File(options.inputfile)
        else:
            raise IOError("Genome filename expected to end with ''.vcf'," +
                          " '.vcf.gz', or '.vcf.bz2'.")
    else:
        sys.stderr.write("Provide input VCF file\n")
        parser.print_help()
        sys.exit(1)

    if options.build and options.build in ['b37', 'b38']:
        build = options.build
    else:
        raise IOError("Input VCF genome build must be 'b37' or 'b38'.")

    if (not (options.clinvarfile or options.clinvardir) or
            (options.clinvarfile and options.clinvardir)):
        sys.stderr.write("Please provide either a ClinVar file or directory.")
        parser.print_help()
        sys.exit(1)
    if options.clinvarfile:
        clinvarfilename = options.clinvarfile
    elif options.clinvardir:
        clinvarfilename = get_latest_vcf_file(target_dir=options.clinvardir,
                                              build=build)
    if clinvarfilename.endswith('.vcf'):
        input_clinvar_file = open(options.clinvarfile)
    elif clinvarfilename.endswith('.vcf.gz'):
        input_clinvar_file = gzip.open(clinvarfilename)
    elif clinvarfilename.endswith('.vcf.bz2'):
        input_clinvar_file = bz2.BZ2File(clinvarfilename)
    else:
        raise IOError("ClinVar filename expected to end with '.vcf'," +
                      " '.vcf.gz', or '.vcf.bz2'.")

    if options.type not in ['csv', 'json']:
        raise IOError("Not a valid report type, must be 'csv' or 'json'.")
    if options.type == "csv":
        csv_report(input_genome_file=input_genome_file,
                   input_clinvar_file=input_clinvar_file,
                   build=build,
                   version=version)
    elif options.type == "json":
        notes_json = {}
        if options.notes:
            notes_json["parameter"] = options.notes
            try:
                notes_json = json.loads(options.notes)
            except:
                sys.stderr.write("Could not parse JSON notes field\n")
        json_report(input_genome_file=input_genome_file,
                    input_clinvar_file=input_clinvar_file,
                    build=build,
                    notes=notes_json,
                    version=version)

if __name__ == "__main__":
    main()
