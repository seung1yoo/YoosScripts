### server : dragon

def main():
    import os
    import glob
    for afa in glob.glob('../00.RawData/*.fa'):
        #print afa.split('/')[-1].split('.')[0]
        sample = afa.split('/')[-1].split('.')[0]
        cmd = 'python /BiO/sgpark/MetagenomePipeline/TaxonomyProfiling.py {0} {1} 3'.format(afa, sample)
        print cmd
        #os.system(cmd)

if __name__=='__main__':
    main()

