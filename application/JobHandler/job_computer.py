#!/hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python
import sys
import gzip
import simplejson as json


class ConfigError(Exception):
    # always exit program if error occurs
    def __init__(self, message):
        super(ConfigError, self).__init__(message)


class ComputeJobParam:

    def __init__(self, path, header_ix, skip, chsize, subchsize, method):
        self.path = path
        self.header = header_ix
        self.skip = skip
        self.chsize = chsize
        self.subchsize = subchsize
        self.method = method
        self.header = self.__header_check()

    def compute(self):
        self.__cnt_lines()
        if self.method == 'single':
            return 1
        n = self.__compute_n(self.chsize, self.lines)
        if self.method == 'subparallel':
            if not self.subchsize:
                raise Exception('Subchunk size should be given in config')
        return n

    def __header_check(self):
        if isinstance(self.header, list):
            return self.skip
        else:
            return None

    def __cnt_lines(self):
        count = 0
        line = True
        with self.__open_up(self.path) as fn:
            while line:
                line = fn.readline()
                if line:
                    count += 1
        head = 1 if isinstance(self.header, int) else 0
        lines = (count - self.skip - head)
        self.lines = lines

    def __compute_n(self, num, denom):
        if num > denom:
            raise ConfigError('Chunk size is larger than number of lines')
        n = int(denom/num)
        if denom % num > 0:
            n += 1
        return n

    def __open_up(self, fn):
        try:
            with gzip.open(fn) as gzfn:
                gzfn.read(2)
            return gzip.open(fn)
        except OSError:
            info = str(sys.exc_info())
            if 'gzipped' in info:
                return open(fn)
            elif 'FileNotFound' in info:
                raise ConfigError('No config file found at path {}'.format(fn))
