__author__ = 'armartin'
from optparse import  OptionParser
from collections import defaultdict, OrderedDict, Callable

USAGE = """
mask_bed.py --bed
            --mask
            --out
"""
parser = OptionParser(USAGE)

parser.add_option('--bed')
parser.add_option('--mask')
parser.add_option('--out')

(options, args) = parser.parse_args()

mask = open(options.mask)
bed = open(options.bed)
out = open(options.out, 'w')

chrs = map(str, range(1,23))
chrs.append('X')

class DefaultOrderedDict(OrderedDict):
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
            not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))
    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                        OrderedDict.__repr__(self))

#masked files need to be ordered so bed can interated through properly
mask_dict = DefaultOrderedDict(OrderedDict)
for line in mask:
    line = line.strip().split()
    mask_dict[line[0]][int(line[1])] = line

last_chr = 0
out = open(options.out, 'w')
while True:
    line = bed.readline().strip().split()
    if not line: break
    try:
        current_chr = int(line[0])
    except ValueError:
        current_chr = 23
    if line[0] in mask_dict:
        for mask_start in mask_dict[line[0]]:
            #mask is the same tract
            if int(line[1]) == mask_start and int(line[2]) == int(mask_dict[line[0]][mask_start][2]):
                out.write('\t'.join(mask_dict[line[0]][mask_start]) + '\n')
            #mask is within a single tract
            elif int(line[1]) <= mask_start and int(line[2]) >= int(mask_dict[line[0]][mask_start][2]):
                if int(line[1]) != mask_start:
                    tract1 = list(line)
                    tract1[2] = str(mask_start-1)
                    tract1[5] = mask_dict[line[0]][mask_start][4]
                    out.write('\t'.join(tract1) + '\n')
                out.write('\t'.join(mask_dict[line[0]][mask_start]) + '\n')
                if int(line[2]) != mask_dict[line[0]][mask_start][5]:
                    tract2 = list(line)
                    tract2[1] = str(int(mask_dict[line[0]][mask_start][2])+1)
                    tract2[4] = mask_dict[line[0]][mask_start][5]
                    tract2[5] = line[5]
                    out.write('\t'.join(tract2) + '\n')
            #mask is within two tracts
            elif int(line[1]) <= mask_start and int(line[2]) <= int(mask_dict[line[0]][mask_start][2]) and int(line[2]) >= mask_start:
                if int(line[1]) != mask_start:
                    tract1 = list(line)
                    tract1[2] = str(mask_start-1)
                    tract1[5] = mask_dict[line[0]][mask_start][4]
                    out.write('\t'.join(tract1) + '\n')
                    #otherwise no tract 1 because mask subsumes it
                out.write('\t'.join(mask_dict[line[0]][mask_start]) + '\n')
                line = bed.readline().strip().split() #get the next tract
                if int(line[2]) != mask_dict[line[0]][mask_start][5]:
                    tract2 = list(line)
                    tract2[1] = str(int(mask_dict[line[0]][mask_start][2])+1)
                    tract2[4] = mask_dict[line[0]][mask_start][5]
                    tract2[5] = line[5]
                    out.write('\t'.join(tract2) + '\n')
                    #otherwise no tract 2 because mask subsumes it
            elif last_chr != current_chr and int(line[1]) > mask_start and int(line[1]) < mask_dict[line[0]][mask_start][2]:
                out.write('\t'.join(mask_dict[line[0]][mask_start]) + '\n')
                if len(mask_dict[line[0]]) > 1:
                    continue
                else:
                    out.write('\t'.join(line) + '\n')
            #elif last_chr != current_chr and int(line[1]) > mask_start and int(line[1]) > mask_dict[line[0]][mask_start][2]:
            #    print 'beginning2 ' + str(mask_start) + ' ' + mask_dict[line[0]][mask_start][2]
            else: #not within mask
                out.write('\t'.join(line) + '\n')
            #tract ends in middle of mask
            last_chr = current_chr
    else:
        out.write('\t'.join(line) + '\n')
out.close()