from __future__ import unicode_literals

from collections import OrderedDict
import json

from .common import Allele, VCFLine


CLNSIG_INDEX = {
    '0': "unknown",
    '1': "untested",
    '2': "non-pathogenic",
    '3': "probably non-pathogenic",
    '4': "probably pathogenic",
    '5': "pathogenic",
    '6': "affecting drug response",
    '7': "affecting histocompatibility",
    '255': "other"}


class ClinVarAllele(Allele):
    """Store ClinVar data relating to one allele."""
    def __init__(self, *args, **kwargs):
        """
        Initialize ClinVarAllele object

        A ClinVarAllele is an allele for a genomic position that has data
        from ClinVar associated with it.

        Required arguments:
        sequence:  String of DNA letters (A, C, G, or T) for the allele;
                   may be empty (to represent a deletion)
        frequency: Preferred allele frequency
        alleleid:  ClinVar Allele ID
        clnhgvs:   HGVS nomenclature for this allele
        clnsig:    ClinVar clinical significance
        clndn:     ClinVar disease name
        clndisdb:  Database IDs of disease database entries (tag-value pairs)
        clnvi:     Database IDs of clinical sources (tag-value pairs)
        """
        (self.clnalleleid, self.hgvs, self.clnsig,
         self.clndn, self.clndisdb, self.clnvi) = [
            kwargs[x] for x in
            ['alleleid', 'clnhgvs', 'clnsig', 'clndn', 'clndisdb', 'clnvi']]
        super(ClinVarAllele, self).__init__(*args, **kwargs)

    def as_dict(self, *args, **kwargs):
        """Return ClinVarAllele data as dict object."""
        self_as_dict = super(ClinVarAllele, self).as_dict(*args, **kwargs)
        self_as_dict['hgvs'] = self.hgvs
        self_as_dict['clnalleleid'] = self.clnalleleid
        self_as_dict['clnsig'] = self.clnsig
        self_as_dict['clndn'] = self.clndn
        self_as_dict['clndisdb'] = self.clndisdb
        self_as_dict['clnvi'] = self.clnvi
        return self_as_dict


class ClinVarVCFLine(VCFLine):
    """Store ClinVar data from a VCF line."""

    def __init__(self, *args, **kwargs):
        """Initialize ClinVarVCFLine with VCF line"""
        kwargs['skip_info'] = False
        super(ClinVarVCFLine, self).__init__(self, *args, **kwargs)

    def as_dict(self):
        """Dict representation of parsed ClinVar VCF line"""
        return {'chrom': self.chrom,
                'start': self.start,
                'ref_allele': self.ref_allele,
                'alt_alleles': self.alt_alleles,
                'info': self.info,
                'alleles': [x.as_dict() for x in self.alleles]}

    def _parse_frequencies(self):
        """Parse frequency data in ClinVar VCF"""
        frequencies = OrderedDict([
            ('EXAC', 'Unknown'),
            ('ESP', 'Unknown'),
            ('TGP', 'Unknown')])
        pref_freq = 'Unknown'
        for source in frequencies.keys():
            freq_key = 'AF_' + source
            if freq_key in self.info:
                frequencies[source] = self.info[freq_key]
                if pref_freq == 'Unknown':
                    pref_freq = frequencies[source]
        return pref_freq, frequencies

    def _parse_allele_data(self):
        """Parse alleles for ClinVar VCF, overrides parent method."""

        # Get allele frequencies if they exist.
        pref_freq, frequencies = self._parse_frequencies()

        info_clnvar_single_tags = ['ALLELEID', 'CLNSIG', 'CLNHGVS']
        cln_data = {x.lower(): self.info[x] if x in self.info else None
                    for x in info_clnvar_single_tags}
        cln_data.update(
            {'clndisdb': [x.split(',') for x in
                          self.info['CLNDISDB'].split('|')]
             if 'CLNDISDB' in self.info else []})
        cln_data.update({'clndn': self.info['CLNDN'].split('|') if
                         'CLNDN' in self.info else []})
        cln_data.update({'clnvi': self.info['CLNVI'].split(',')
                        if 'CLNVI' in self.info else []})

        try:
            sequence = self.alt_alleles[0]
        except IndexError:
            sequence = self.ref_allele

        allele = ClinVarAllele(frequency=pref_freq, sequence=sequence,
                               **cln_data)

        # A few ClinVar variants are only reported as a combination with
        # other variants, and no single-variant effect is proposed. Skip these.
        if not cln_data['clnsig']:
            return []

        return [allele]
