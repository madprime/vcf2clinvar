from __future__ import unicode_literals
import json
import re


CHROM_INDEX = {
    '1': 1, '2': 2, '3': 3, '4': 4, '5': 5,
    '6': 6, '7': 7, '8': 8, '9': 9, '10': 10,
    '11': 11, '12': 12, '13': 13, '14': 14, '15': 15,
    '16': 16, '17': 17, '18': 18, '19': 19, '20': 20,
    '21': 21, '22': 22, 'X': 23, 'Y': 24, 'M': 25, 'MT': 25,
    'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5,
    'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10,
    'chr11': 11, 'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15,
    'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 'chr20': 20,
    'chr21': 21, 'chr22': 22, 'chrX': 23, 'chrY': 24, 'chrM': 25,
}

REV_CHROM_INDEX = {
    1: 'chr1', 2: 'chr2', 3: 'chr3', 4: 'chr4', 5: 'chr5',
    6: 'chr6', 7: 'chr7', 8: 'chr8', 9: 'chr9', 10: 'chr10',
    11: 'chr11', 12: 'chr12', 13: 'chr13', 14: 'chr14', 15: 'chr15',
    16: 'chr16', 17: 'chr17', 18: 'chr18', 19: 'chr19', 20: 'chr20',
    21: 'chr21', 22: 'chr22', 23: 'chrX', 24: 'chrY', 25: 'chrM',
}


class Allele(object):
    """Store data relating to one allele.

    Instance attributes:
    sequence:  (string) The DNA sequence of this allele, may only contain the
               characters A, C, G, T, or N.
    frequency: (string) String representation of a number between 0 and 1,
               or "Unknown".
    """
    def __init__(self, *args, **kwargs):
        """
        Initialize Allele object

        Required arguments:
        sequence:  Short string of DNA letters (ACGT) for the allele.
                   May be empty (to represent a deletion).

        Optional arguments:
        frequency: a string representation of a float between 0 and 1
        """
        sequence = kwargs['sequence']
        if 'frequency' in kwargs:
            frequency = kwargs['frequency']
        else:
            frequency = 'Unknown'

        if not (re.match(r'^[ACGTN]*$', sequence) or
                re.match(r'^<.*>$', sequence)):
            raise ValueError("Allele sequence isn't a standard DNA sequence")
        self.sequence = sequence
        if frequency:
            try:
                if (float(frequency) < 0.0 or
                        float(frequency) > 1.0):
                    raise ValueError('Allele frequency not between 0 and 1')
            except ValueError:
                if not frequency == 'Unknown':
                    raise ValueError('Allele frequency must be a number ' +
                                     'between 0 and 1, or the string ' +
                                     '"Unknown".')
            self.frequency = frequency

    def __unicode__(self):
        """Print Allele object as dict object data."""
        return self.as_json()

    def __str__(self):
        """Print Allele object as dict object data."""
        return self.as_json()

    def as_dict(self):
        """Return Allele data as dict object."""
        self_as_dict = dict()
        self_as_dict['sequence'] = self.sequence
        if hasattr(self, 'frequency'):
            self_as_dict['frequency'] = self.frequency
        return self_as_dict

    def as_json(self):
        """Print Allele object as JSON."""
        return json.dumps(self.as_dict())


class VCFLine(object):
    """Process data from a VCF line."""

    def __init__(self, *args, **kwargs):
        """Store data from a VCF line."""
        vcf_line = kwargs['vcf_line']
        skip_info = ('skip_info' in kwargs and kwargs['skip_info'])

        vcf_fields = vcf_line.strip().split('\t')
        self.chrom = vcf_fields[0]
        self.start = int(vcf_fields[1])
        self.ref_allele = vcf_fields[3]
        if vcf_fields[4] == '.':
            self.alt_alleles = []
        else:
            self.alt_alleles = vcf_fields[4].split(',')
        if not skip_info:
            self.info = self._parse_info(vcf_fields[7])
        self.alleles = self._parse_allele_data()

    def _parse_allele_data(self):
        """Create list of Alleles from VCF line data"""
        return [Allele(sequence=x) for x in
                [self.ref_allele] + self.alt_alleles]

    def _parse_info(self, info_field):
        """Parse the VCF info field"""
        info = dict()
        for item in info_field.split(';'):
            # Info fields may be "foo=bar" or just "foo".
            # For the first case, store key "foo" with value "bar"
            # For the second case, store key "foo" with value True.
            info_item_data = item.split('=')
            # If length is one, just store as a key with value = true.
            if len(info_item_data) == 1:
                info[info_item_data[0]] = True
            elif len(info_item_data) == 2:
                info[info_item_data[0]] = info_item_data[1]
        return info

    def __str__(self):
        """String representation of parsed VCF data"""
        return json.dumps(self.as_dict(), ensure_ascii=True)

    def as_dict(self):
        """Dict representation of parsed VCF data"""
        self_as_dict = {'chrom': self.chrom,
                        'start': self.start,
                        'ref_allele': self.ref_allele,
                        'alt_alleles': self.alt_alleles,
                        'alleles': [x.as_dict() for x in self.alleles]}
        try:
            self_as_dict['info'] = self.info
        except AttributeError:
            pass
        return self_as_dict

    @staticmethod
    def get_pos(vcf_line):
        """
        Very lightweight parsing of a vcf line to get position.

        Returns a dict containing:
        'chrom': index of chromosome (int), indicates sort order
        'pos': position on chromosome (int)
        """
        if not vcf_line:
            return None
        vcf_data = vcf_line.strip().split('\t')
        return_data = dict()
        return_data['chrom'] = CHROM_INDEX[vcf_data[0]]
        return_data['pos'] = int(vcf_data[1])
        return return_data
