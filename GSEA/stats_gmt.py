
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def load_gmt(gmt):
    unit_count_s = list()
    for line in open(gmt):
        items = line.rstrip('\n').split('\t')
        name = items[0]
        acc = items[1]
        units = items[2:]
        unit_count_s.append(len(units))
    return unit_count_s


def main(args):
    unit_count_s = load_gmt(args.gmt)
    array = np.array(unit_count_s)
    print(array)
    print('mean   : {0}'.format(np.mean(array)))
    print('var    : {0}'.format(np.var(array)))
    print('std    : {0}'.format(np.std(array)))
    print('max    : {0}'.format(np.max(array)))
    print('min    : {0}'.format(np.min(array)))
    print('median : {0}'.format(np.median(array)))
    #
    n, bins, patches = plt.hist(array, bins=100)
    #print(n)
    #print(bins)
    #print(patches)
    plt.savefig('stats_gmt.hist.png')



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('gmt')
    args = parser.parse_args()
    main(args)
