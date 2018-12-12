
import csv
from hgsc_vcf.metainfo import *
from collections import *

class Reader(object):
    def __init__(self, fobj):
        self.fobj = fobj
        self.header = VCFHeader()
        self.header.load(self.fobj)
        self._next = None
    
    @staticmethod
    def parse_info_field(info):
        infos = info.split(';')
        result = OrderedDict()
        for i in infos:
            if '=' in i:
                k, v = i.split('=',1)
                result[k] = v.split(',')
            else:
                result[i] = True # True indicates that the flag is active
        return result

    def peek(self):
        return self._next

    def take(self):
        old = self._next
        try:
            self._next = self.next()
        except StopIteration:
            pass # swallow the error
        return old

    @staticmethod
    def parse_sample(format_keys, slist):
        return OrderedDict(zip(format_keys, [i.split(',') for i in slist]))
    
    def __iter__(self):
        return self

    def next(self):
        line = [c.strip() for c in self.fobj.readline().split('\t')]
        if len(line) < 1 or line[0] == '':
            self._next = None
            raise StopIteration
        try:
            record = OrderedDict()
            for k, v in (
                        ('CHROM', line[0]),
                        ('POS', int(line[1])),
                        ('ID', line[2].split(';')),
                        ('REF', line[3]),
                        ('ALT', line[4].split(',')),
                        ('QUAL', float(line[5]) if line[5] != '.' else '.'),
                        ('FILTER', line[6].split(';')),
                        ('INFO', Reader.parse_info_field(line[7]))
                        ):
                record[k] = v
            if len(line) > 8:
                record['FORMAT'] = line[8].split(':')
                if len(line) > 9:
                    record['SAMPLES'] = OrderedDict(zip(
                        self.header.samples, 
                        [Reader.parse_sample(record['FORMAT'], s.split(':')) for s in line[9:]]
                        ))
            self._next = record
            return record
        except:
            print line
            raise

class Writer(object):
    def __init__(self, fobj, header):
        assert isinstance(header, VCFHeader), "header must be a VCFHeader"
        self.header = header
        self.header_written = False
        self.fobj = fobj

    def write_header(self):
        if self.header_written:
            raise ValueError("Can't write the header twice")
        for h in self.header.headers:
            self.fobj.write(str(h) + '\n')
        header_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        if len(self.header.samples) > 0:
            header_cols.append('FORMAT')
            for s in self.header.samples:
                header_cols.append(s)
        self.fobj.write('#' + '\t'.join(header_cols) + '\n')
        self.header_written = True

    def write_record(self, record):
        if not self.header_written:
            raise ValueError("Must write the header first")
        field_parts = []
        for k, joiner in (('CHROM', None), ('POS', None), ('ID', ';'), ('REF', None), ('ALT', ','), ('QUAL', None), ('FILTER', ';')):
            if joiner:
                try:
                    field_parts.append(joiner.join(record[k]))
                except:
                    print k, joiner, record[k]
                    raise
            else:
                field_parts.append(str(record[k]))
        # info is a bit trickier
        info_parts = []
        for k, v in record['INFO'].items():
            if k == '.' and len(record['INFO']) > 1:
                continue # this is a leftover empty marker
            if isinstance(v, list):
                info_parts.append('%s=%s' % (k, ','.join(v)))
            else:
                info_parts.append(k)
        field_parts.append(';'.join(info_parts))
        if len(self.header.samples) > 0:
            field_parts.append(':'.join(record['FORMAT']))
            for s in self.header.samples:
                sinfo = record['SAMPLES'][s]
                # sinfo is a dict (OrderedDict ideally)
                try:
                    field_parts.append(':'.join([','.join(sinfo[k]) for k in record['FORMAT']]))
                except:
                    print sinfo
                    raise
        self.fobj.write('\t'.join(field_parts) + '\n')


