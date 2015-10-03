from __future__ import unicode_literals
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


class ClinVarRecord(object):
    """
    ClinVar data relating to one record.

    Instance attributes:
    dsdb: ClinVar Record disease database links (array of tuples of strings)
          [("disease-dasabase", "disease-database-id"), ...]
    acc:  ClinVar record accession (string)
    dbn:  ClinVar record disease name (string)
    sig:  ClinVar record significance (string)
    """
    def __init__(self, clndsdb, clndsdbid, clnacc, clndbn, clnsig):
        """
        Initialize ClinVarRecord

        A ClinVar Record is an assertion of effect made by a particular
        submitter. Any given variant may have more than one record, and
        the associated assertions may assess different significance, and may
        be related to different diseases.

        Required arguments:
        clndsdbs:   (string) Disease databases, joined by colon (':')
        clndsdbids: (string) Disease database IDs, joined by colon (':')
        clnacc:     (string) ClinVarRecord accession
        clndbn:     (string) ClinVarRecord disease name
        clnsig:     (string) String representation of a digit; corresponds to
                    keys in CLNSIG_INDEX.
        """
        clndsdbs = clndsdb.split(':')
        clndsdbids = clndsdbid.split(':')
        self.dsdb = [(clndsdbs[i], clndsdbids[i]) for i
                     in range(len(clndsdbs))]
        self.acc = clnacc
        self.dbn = clndbn
        self.sig = clnsig

    def __unicode__(self):
        """Return ClinVarRecord data as JSON-formatted string."""
        return json.dumps(self.as_dict(), ensure_ascii=True)

    def as_dict(self):
        """Return ClinVarRecord data as dict object."""
        return {'dsdb': self.dsdb,
                'acc': self.acc,
                'dbn': self.dbn,
                'sig': self.sig}


class ClinVarAllele(Allele):
    """Store ClinVar data relating to one allele."""
    def __init__(self, *args, **kwargs):
        """
        Initialize ClinVarAllele object

        A ClinVarAllele is an allele for a genomic position that has data
        from ClinVar associated with it. ClinVar data in the VCF appears to
        exist at two levels: the allele, and records. A ClinVar "record"
        describes a reported effect, and any allele may have multiple
        records associated with it. An allele can also have other data
        returned by ClinVar (in the ClinVar VCF this appears to be separate
        to the Records).

        Required arguments:
        sequence:  String of DNA letters (A, C, G, or T) for the allele;
                   may be empty (to represent a deletion).
        records:   list of ClinVarRecord objects associated with this allele
        hgvs:      HGVS nomenclature for this allele
        clnsrcs:   list of ClinVar sources
        clnsrcids: list of IDs for the ClinVar sources

        Optional arguments:
        frequency: a float between 0 and 1, or string saying "Unknown"
        """
        clnsrcs, clnsrcids, clnhgvs, records = [kwargs[x] for x in
                                                ['clnsrcs', 'clnsrcids',
                                                 'clnhgvs', 'records']]
        self.src = [(clnsrcs[i], clnsrcids[i]) for i in range(len(clnsrcs))]
        self.hgvs = clnhgvs
        self.records = records
        super(ClinVarAllele, self).__init__(*args, **kwargs)

    def as_dict(self, *args, **kwargs):
        """Return ClinVarAllele data as dict object."""
        self_as_dict = super(ClinVarAllele, self).as_dict(*args, **kwargs)
        self_as_dict['hgvs'] = self.hgvs
        self_as_dict['src'] = self.src
        self_as_dict['records'] = [x.as_dict() for x in self.records]
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
        given_freqs = self.info['CAF'].rstrip(']').lstrip('[').split(',')
        parsed_freqs = ['Unknown' if x == '.' else x for x in given_freqs]
        return parsed_freqs

    def _parse_clinvar_allele(self, *args, **kwargs):
        """Parse ClinVar records for each allele"""
        cln_data, cln_idx, allele_idx = [kwargs[x] for x in
                                         ['cln_data', 'cln_idx', 'allele_idx']]
        if 'frequency' in kwargs:
            frequency = kwargs['frequency']
        else:
            frequency = 'Unknown'

        if allele_idx == 0:
            sequence = self.ref_allele
        else:
            sequence = self.alt_alleles[allele_idx - 1]
        clnsrcs = cln_data['CLNSRC'][cln_idx]
        clnsrcids = cln_data['CLNSRCID'][cln_idx]
        clnhgvs = cln_data['CLNHGVS'][cln_idx][0]

        # Process all the ClinVar records for this allele.
        records = []
        for record_idx in range(len(cln_data['CLNACC'][cln_idx])):
            try:
                record = ClinVarRecord(
                    clndsdb=cln_data['CLNDSDB'][cln_idx][record_idx],
                    clndsdbid=cln_data['CLNDSDBID'][cln_idx][record_idx],
                    clnacc=cln_data['CLNACC'][cln_idx][record_idx],
                    clndbn=cln_data['CLNDBN'][cln_idx][record_idx],
                    clnsig=cln_data['CLNSIG'][cln_idx][record_idx])
            except IndexError:
                # Skip inconsintent entries. At least one line in the
                # ClinVar VCF as of 2014/06 has inconsistent CLNSIG and
                # CLNACC information (rs799916).
                return self._parse_allele(*args, **kwargs)
            records.append(record)
        return ClinVarAllele(sequence=sequence,
                             clnhgvs=clnhgvs,
                             clnsrcs=clnsrcs,
                             clnsrcids=clnsrcids,
                             records=records,
                             frequency=frequency)

    def _parse_allele(self, *args, **kwargs):
        """Create an Allele, with optional frequency data."""
        allele_idx = kwargs['allele_idx']
        try:
            frequency = kwargs['frequency']
        except KeyError:
            frequency = 'Unknown'

        if allele_idx == 0:
            sequence = self.ref_allele
        else:
            sequence = self.alt_alleles[allele_idx - 1]
        if frequency:
            return Allele(sequence=sequence,
                          frequency=frequency)
        else:
            return Allele(sequence=sequence)

    def _parse_allele_data(self):
        """Parse alleles, overrides parent method."""
        # Get allele frequencies if they exist.
        frequencies = []
        if 'CAF' in self.info:
            frequencies = self._parse_frequencies()

        # CLNALLE describes which allele ClinVar data correspond to.
        clnalle_keys = [int(x) for x in self.info['CLNALLE'].split(',')]
        info_clnvar_tags = ['CLNDSDB', 'CLNDSDBID', 'CLNACC', 'CLNDBN',
                            'CLNSIG', 'CLNHGVS', 'CLNSRC', 'CLNSRCID']
        # Clinvar data is split first by comma, then by pipe.
        cln_data = {x: [y.split('|') for y in self.info[x].split(',')]
                    for x in info_clnvar_tags if x}

        # Iterate over all alleles, if index is in clnallele_keys then
        # create a ClinVarAllele, otherwise create an Allele.
        alleles = []
        for i in range(len(self.alt_alleles) + 1):
            if i in clnalle_keys and frequencies:
                cln_idx = clnalle_keys.index(i)
                allele = self._parse_clinvar_allele(allele_idx=i,
                                                    cln_idx=cln_idx,
                                                    cln_data=cln_data,
                                                    frequency=frequencies[i])
            elif i in clnalle_keys:
                cln_idx = clnalle_keys.index(i)
                allele = self._parse_clinvar_allele(allele_idx=i,
                                                    cln_idx=cln_idx,
                                                    cln_data=cln_data)
            elif frequencies:
                allele = self._parse_allele(allele_idx=i,
                                            frequency=frequencies[i])
            else:
                allele = self._parse_allele(allele_idx=i)
            alleles.append(allele)
        return alleles
