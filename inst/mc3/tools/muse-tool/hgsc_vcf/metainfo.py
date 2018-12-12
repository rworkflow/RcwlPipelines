
import re
from collections import *

class HeaderLine(object):
    SIMPLE = 1
    COMPLEX = 2
    _RESERVED_KEYS = ['ID', 'Number', 'Type']
    def __init__(self, key, type, fields):
        assert isinstance(key, basestring), "key must be a string type"
        assert isinstance(type, int), "Type must be one of HeaderLine.SIMPLE or HeaderLine.COMPLEX"
        assert isinstance(fields, OrderedDict), "fields must be a mapping of header line info data"
        self.type = type
        self.fields = fields
        self.key = key

    def __getitem__(self, key):
        return self.fields[key]
    
    def __str__(self):
        if self.type == HeaderLine.COMPLEX:
            return '##%(key)s=<%(fields)s>' % {
                    'key': self.key,
                    'fields': ','.join([HeaderLine.format_field(k, v) for k, v in self.fields.items()])
                    }
        else:
            return '##%(key)s=%(fields)s' % {
                    'key': self.key,
                    'fields': ','.join([str(v) for v in self.fields.values()])
                    }
    
    @staticmethod
    def format_field(k, v):
        if k in HeaderLine._RESERVED_KEYS:
            return '%s=%s' % (k, str(v))
        else:
            return '%s="%s"' % (k, str(v))

class ComplexHeaderLine(HeaderLine):
    BQUOTE = '__BQUOTE__'
    PATTERN = re.compile(r'''((?:[^,"']|"[^"]*"|'[^']*')+)''')
    COMPLEX = re.compile(r'''##(.*?)=<(.*)>$''')
    SIMPLE = re.compile(r'''##(.*?)=(.*)$''')
    QUOTED = re.compile(r'''["'](.*?)["']''')
    def __init__(self, key = None, fields = None, string = None, hltype = HeaderLine.COMPLEX):
        if fields is not None:
            assert isinstance(key, basestring), "Must specify a key with fields option"
            assert isinstance(fields, OrderedDict), "Fields must either be an OrderedDict or None"
            assert string is None, "You may not specify a string for parsing and fields at the same time"
            HeaderLine.__init__(self, key, hltype, fields)
        elif string is not None:
            assert isinstance(string, basestring), "string must be a string type"
            key, fields, hltype = ComplexHeaderLine.parse_string(string)
            HeaderLine.__init__(self, key, hltype, fields)
        else:
            raise ValueError("Can't create a header line with no fields without explicitly specifying string")
    
    @staticmethod
    def parse_string(string):
        string = string.replace('\\"', ComplexHeaderLine.BQUOTE).strip()
        l_m = re.match(ComplexHeaderLine.COMPLEX, string)
        if l_m:
            k, f_s = l_m.group(1,2)
            items = []
            for f_i in ComplexHeaderLine.PATTERN.split(str(f_s))[1::2]: # start with index = 1 and go every other one
                if '=' not in f_i:
                    f_i = f_i + '=True'
                k_i, v_i = f_i.split('=', 1)
                v_iq = re.match(ComplexHeaderLine.QUOTED, v_i)
                if v_iq:
                    v_i = v_iq.group(1)
                items.append((k_i, v_i.replace(ComplexHeaderLine.BQUOTE, '\\"')))
            f = OrderedDict(items)
            return k, f, HeaderLine.COMPLEX
        l_m = re.match(ComplexHeaderLine.SIMPLE, string)
        if l_m:
            k, f_s = l_m.group(1,2)
            return k, OrderedDict(value = f_s.replace(ComplexHeaderLine.BQUOTE, '\\"')), HeaderLine.SIMPLE


class VCFHeader(object):
    def __init__(self):
        self.headers = []
        self.samples = []
    @staticmethod
    def _header_line_matches(h, htype, id):
        if h.key != htype:
            return False
        if not id:
            return True
        return h.fields['ID'] == id


    def get_headers(self, htype, id = None):
        for h in self.headers:
            if h.key == htype:
                if id:
                    if h.fields.get('ID', None) == id:
                        yield h
                else:
                    yield h

    def add_header(self, string):
        h = ComplexHeaderLine(string = string)
        self.headers.append(h)

    def remove_header(self, htype, id = None):
        new_headers = [h for h in self.headers if not VCFHeader._header_line_matches(h, htype, id)]
        self.headers = new_headers

    def set_headers(self, new_headers):
        assert isinstance(new_headers, list), "headers must be a list"
        for h in new_headers:
            assert isinstance(h, HeaderLine), "each header must be a HeaderLine"
        self.headers = new_headers

    def get_format_keys(self):
        format_header_lines = self.get_headers('FORMAT')
        return [f.fields['ID'] for f in format_header_lines]

    def get_info_keys(self):
        return [f.fields['ID'] for f in self.get_headers('INFO')]

    def load(self, fobj):
        line = fobj.readline().strip()
        if not line.startswith('#'):
            raise ValueError("The first line of any VCF should begin with a #")
        while not line.startswith('#CHROM'):
            self.add_header(line)
            line = fobj.readline().strip()
        header_cols = [c.strip() for c in line.replace('#', '').split('\t')]
        if len(header_cols) > 9:
            self.samples = header_cols[9:]
