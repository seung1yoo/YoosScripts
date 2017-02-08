

def main(args):
    out = open(args.outxls, 'w')
    for line in open(args.inxls):
        new_items = []
        items = line.rstrip('\n').split('\t')
        #print(items[:11])
        new_items.extend(items[:11])
        tgItems = items[11:]
        #print(tgItems)
        tgItems = ['0.00' if item in ['-'] else item for item in tgItems]
        #print(tgItems)
        new_items.extend(tgItems)
        out.write('{0}\n'.format('\t'.join(new_items)))
    out.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inxls')
    parser.add_argument('-o', '--outxls')
    args = parser.parse_args()
    main(args)
