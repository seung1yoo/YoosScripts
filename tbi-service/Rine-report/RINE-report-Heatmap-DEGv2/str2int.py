

def main():
    out = open('genes.90.addDEG.pathogenic.int.xls', 'w')
    for line in open('genes.90.addDEG.pathogenic.xls'):
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
    main()
