#!/usr/bin/python
#
"""Tools for parsing and matching VCF files"""

from .clinvar import CLNSIG_INDEX, ClinVarVCFLine
from .common import VCFLine
from .genome import GenomeVCFLine


def _next_line(filebuffer):
    try:
        next_line = filebuffer.readline()
    except AttributeError:
        next_line = filebuffer.next()
    try:
        next_line = next_line.decode('utf-8')
        return next_line
    except AttributeError:
        return next_line


def match_to_clinvar(genome_file, clin_file):
    """
    Match a genome VCF to variants in the ClinVar VCF file

    Acts as a generator, yielding tuples of:
    (ClinVarVCFLine, ClinVarAllele, zygosity)

    'zygosity' is a string and corresponds to the genome's zygosity for that
    ClinVarAllele. It can be either: 'Het' (heterozygous), 'Hom' (homozygous),
    or 'Hem' (hemizygous, e.g. X chromosome in XY individuals).
    """
    clin_curr_line = _next_line(clin_file)
    genome_curr_line = _next_line(genome_file)

    # Ignores all the lines that start with a hashtag
    while clin_curr_line.startswith('#'):
        clin_curr_line = _next_line(clin_file)
    while genome_curr_line.startswith('#'):
        genome_curr_line = _next_line(genome_file)

    # Advance through both files simultaneously to find matches
    while clin_curr_line and genome_curr_line:

        # Advance a file when positions aren't equal.
        clin_curr_pos = VCFLine.get_pos(clin_curr_line)
        genome_curr_pos = VCFLine.get_pos(genome_curr_line)
        try:
            if clin_curr_pos['chrom'] > genome_curr_pos['chrom']:
                genome_curr_line = _next_line(genome_file)
                continue
            elif clin_curr_pos['chrom'] < genome_curr_pos['chrom']:
                clin_curr_line = _next_line(clin_file)
                continue
            if clin_curr_pos['pos'] > genome_curr_pos['pos']:
                genome_curr_line = _next_line(genome_file)
                continue
            elif clin_curr_pos['pos'] < genome_curr_pos['pos']:
                clin_curr_line = _next_line(clin_file)
                continue
        except StopIteration:
            break

        # If we get here, start positions match.
        # Look for allele matching.
        genome_vcf_line = GenomeVCFLine(vcf_line=genome_curr_line,
                                        skip_info=True)
        # We can skip if genome has no allele information for this point.
        if not genome_vcf_line.genotype_allele_indexes:
            genome_curr_line = _next_line(genome_file)
            continue

        # Match only if ClinVar and Genome ref_alleles match.
        clinvar_vcf_line = ClinVarVCFLine(vcf_line=clin_curr_line)
        if not genome_vcf_line.ref_allele == clinvar_vcf_line.ref_allele:
            try:
                genome_curr_line = _next_line(genome_file)
                clin_curr_line = _next_line(clin_file)
                continue
            except StopIteration:
                break

        # Determine genome alleles and zygosity. Zygosity is assumed to be one
        # of: heterozygous, homozygous, or hemizygous.
        genotype_allele_indexes = genome_vcf_line.genotype_allele_indexes
        genome_alleles = [genome_vcf_line.alleles[x] for
                          x in genotype_allele_indexes]
        if len(genome_alleles) == 1:
            zygosity = 'Hem'
        elif len(genome_alleles) == 2:
            if genome_alleles[0].sequence == genome_alleles[1].sequence:
                zygosity = 'Hom'
                genome_alleles = [genome_alleles[0]]
            else:
                zygosity = 'Het'
        else:
            raise ValueError('This code only expects to work on genomes ' +
                             'with one or two alleles called at each ' +
                             'location. The following line violates this:' +
                             str(genome_vcf_line))

        # Look for matches to ClinVar alleles.
        for genome_allele in genome_alleles:
            for allele in clinvar_vcf_line.alleles:
                if genome_allele.sequence == allele.sequence:
                    # The 'records' attribute is specific to ClinVarAlleles.
                    if hasattr(allele, 'records'):
                        yield (genome_vcf_line, allele, zygosity)

        # Done matching, move on.
        try:
            genome_curr_line = _next_line(genome_file)
            clin_curr_line = _next_line(clin_file)
        except StopIteration:
            break
