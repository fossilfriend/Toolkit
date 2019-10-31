class FakeSecHead(object):
    '''
    puts a fake section into the gus.config file so that ConfigParser can
    be used to retrieve db access information

    workaround for the INI-style headings required by ConfigParser
    see https://stackoverflow.com/questions/2819696/parsing-properties-file-in-python
    '''
    def __init__(self, fp):
        self.fp = fp
        self.sechead = '[section]\n'

    def readline(self):
        if self.sechead:
            try: 
                return self.sechead
            finally: 
                self.sechead = None
        else: 
            return self.fp.readline()
