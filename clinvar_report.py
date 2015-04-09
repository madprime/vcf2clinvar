#!/usr/bin/python
#

from __future__ import unicode_literals
from argparse import ArgumentParser
import bz2
import csv
import gzip
import json
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
import sys

from vcf2clinvar import match_to_clinvar
from vcf2clinvar.common import REV_CHROM_INDEX

# Don't stack dump on keyboard ctrl-c or on premature
# termination of output stream (say from piping output
# through head).
#
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)

def main():
    """
    Parse command line argument and
    output appropriate file type (csv or JSON)
    """
    parser = ArgumentParser()

    parser.add_argument(
        "-C", "--clinvar", dest="clinvar",
        help="ClinVar VCF file", metavar="CLINVAR")
    parser.add_argument(
        "-i", "--input", dest="inputfile",
        help="Input VCF file", metavar="INPUT")
    parser.add_argument(
        "-F", "--output-format", dest="format",
        help="Output format (currently 'csv' or 'json')",
        metavar="FORMAT")
    parser.add_argument(
        "-V", "--schema-version", dest="schema_version",
        help="Version to include report (JSON only)",
        metavar="OUTVERSION")
    parser.add_argument(
        "-n", "--notes", dest="notes",
        help="Notes, as a JSON string, to include in report (JSON only)",
        metavar="NOTES")
    parser.add_argument(
        "-g", "--genome-build", dest="build",
        help="Genome build to include in report (JSON only)",
        metavar="GENOMEBUILD")
    options = parser.parse_args()

    if sys.stdin.isatty():
        if options.inputfile:
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
    else:
        input_genome_file = sys.stdin

    if options.clinvar:
        if options.clinvar.endswith('.vcf'):
            input_clinvar_file = open(options.clinvar)
        elif options.clinvar.endswith('.vcf.gz'):
            input_clinvar_file = gzip.open(options.clinvar)
        elif options.clinvar.endswith('.vcf.bz2'):
            input_clinvar_file = bz2.BZ2File(options.clinvar)
        else:
            raise IOError("ClinVar filename expected to end with '.vcf'," +
                          " '.vcf.gz', or '.vcf.bz2'.")
    else:
        sys.stderr.write("Provide ClinVar VCF file\n")
        parser.print_help()
        sys.exit(1)

    output_format = "csv"
    if options.format:
        if options.format == "csv":
            output_format = "csv"
        elif options.format == "json":
            output_format = "json"

    if output_format == "csv":
        csv_out = csv.writer(sys.stdout)
        header = ("Chromosome", "Position", "Name", "Significance",
                  "Frequency", "Zygosity", "ACC URL")
        csv_out.writerow(header)

    metadata = {}
    metadata["notes"] = options.clinvar

    build = "unknown"
    if options.build:
        build = options.build
    metadata["genome_build"] = build

    notes_json = {}
    if options.notes:
        notes_json["parameter"] = options.notes
        try:
            notes_json = json.loads(options.notes)
        except:
            sys.stderr.write("Could not parse JSON notes field\n")

    json_report = {}
    json_report["schema_version"] = options.schema_version
    json_report["notes"] = notes_json
    json_report["metadata"] = metadata
    json_report["variants"] = []

    matching = match_to_clinvar(input_genome_file, input_clinvar_file)
    for var in matching:
        chrom = var[0]
        pos = var[1]
        ref_allele = var[2]
        alt_allele = var[3]
        name_acc = var[4]
        allele_freq = var[5]
        zygosity = var[6]

        for spec in name_acc:
            ele = {}
            ele["chrom"] = REV_CHROM_INDEX[chrom]
            ele["pos"] = pos
            ele["ref_allele"] = ref_allele
            ele["alt_allele"] = alt_allele
            ele["allele_freq"] = allele_freq
            ele["zygosity"] = zygosity

            url = 'http://www.ncbi.nlm.nih.gov/clinvar/' + spec[0]
            name = spec[1]
            clnsig = spec[2]

            ele["acc_url"] = url
            ele["name"] = name
            ele["clinical_significance"] = clnsig

            if output_format == 'json':
                json_report["variants"].append(ele)

            if output_format == "csv":
                data = (chrom, pos, name, clnsig, allele_freq, zygosity, url)
                csv_out.writerow(data)

    if output_format == "json":
        print(json.dumps(json_report))


if __name__ == "__main__":
    main()
