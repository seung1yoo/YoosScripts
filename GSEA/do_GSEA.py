
def designParser(design):
    pair = design.replace(',','_versus_')
    return pair

def command_maker(jar, chip, cls, exp, gmt, outdir, rpt_label, pair, metric):
    import os
    if not os.path.isdir(outdir):
        os.system('mkdir -p {0}'.format(outdir))
    #
    command = '''java -Xmx2G -cp {0} xtools.gsea.Gsea \
-res {1} \
-cls {2}#{3} \
-gmx {4} \
-chip {5} \
-out {6} \
-rpt_label {7} \
-metric {8} \
-collapse false \
-mode Max_probe \
-norm meandiv \
-nperm 1000 \
-permute gene_set \
-rnd_type no_balance \
-scoring_scheme weighted \
-sort real \
-order descending \
-include_only_symbols true \
-make_sets true \
-median false \
-num 100 \
-plot_top_x 100 \
-rnd_seed timestamp \
-save_rnd_lists false \
-set_max 500 \
-set_min 10 \
-zip_report false \
-gui false'''.format(jar, exp, cls, pair, gmt, chip, outdir, rpt_label, metric)
    return command

def execute(command):
    print command
    os.system(command)

def main(args):
    pair = designParser(args.design)
    command = command_maker(args.jar, args.chip, args.cls, args.exp, args.gmt, args.outdir, args.rpt_label, pair, args.metric)
    execute(command)

if __name__=='__main__':
    import os
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--jar')
    parser.add_argument('-ch', '--chip')
    parser.add_argument('-cl', '--cls')
    parser.add_argument('-e', '--exp')
    parser.add_argument('-g', '--gmt')
    parser.add_argument('-o', '--outdir', help='./')
    parser.add_argument('-rl', '--rpt-label', help='DEG001')
    parser.add_argument('-d', '--design', help='control,case')
    parser.add_argument('-m', '--metric', choices=('Diff_of_Classes', 'Signal2Noise'))
    args = parser.parse_args()
    main(args)
