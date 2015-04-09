from __future__ import unicode_literals
import re

from .common import VCFLine


class GenomeVCFLine(VCFLine):
    """Store Genome data from a VCF line."""
    def __init__(self, *args, **kwargs):
        super(GenomeVCFLine, self).__init__(self, *args, **kwargs)
        vcf_line = kwargs['vcf_line']
        vcf_fields = vcf_line.strip().split('\t')
        self.genotype_allele_indexes = self._parse_genotype(vcf_fields)

    def _parse_genotype(self, vcf_fields):
        """Parse genotype from VCF line data"""
        format_col = vcf_fields[8].split(':')
        genome_data = vcf_fields[9].split(':')
        try:
            gt_idx = format_col.index('GT')
        except ValueError:
            return []
        return [int(x) for x in re.split(r'[\|/]', genome_data[gt_idx]) if
                x != '.']
