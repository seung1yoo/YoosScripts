
def designParser(design):
    pair = design.replace(',','_versus_')
    return pair

def command_maker(jar, chip, cls, exp, gmt, outdir, rpt_label, pair, metric, permute):
    import os
    if not os.path.isdir(outdir):
        os.system('mkdir -p {0}'.format(outdir))
    # with jar
    _command = '''java -Xmx1024m -cp {0} xtools.gsea.Gsea \
    -res  {1} \
    -cls  {2}#{3} \
    -gmx  {4} \
    -chip {5} \
    -out  {6} \
    -rpt_label {7} \
    -metric {8} \
    -permute {9} \
    -collapse false -mode Max_probe -norm meandiv -nperm 1000 \
    -rnd_type no_balance -scoring_scheme weighted -sort real \
    -order descending -create_gcts true -create_svgs true \
    -include_only_symbols true -make_sets true -median false \
    -num 100 -plot_top_x 50 -rnd_seed timestamp -save_rnd_lists false \
    -set_max 100 -set_min 5 -zip_report false \
    -gui false'''.format(jar, exp, cls, pair, gmt, chip, outdir, rpt_label, metric, permute)
    #~/YoosScripts/GSEA/GSEA_4.0.2/gsea-cli.sh
    command = '''{0} GSEA \
    -res  {1} \
    -cls  {2}#{3} \
    -gmx  {4} \
    -chip {5} \
    -out  {6} \
    -rpt_label {7} \
    -metric {8} \
    -permute {9} \
    -collapse false -mode Max_probe -norm meandiv -nperm 1000 \
    -rnd_type no_balance -scoring_scheme weighted -sort real \
    -order descending -create_gcts true -create_svgs true \
    -include_only_symbols true -make_sets true -median false \
    -num 100 -plot_top_x 50 -rnd_seed timestamp -save_rnd_lists false \
    -set_max 100 -set_min 5 -zip_report false \
    '''.format(jar, exp, cls, pair, gmt, chip, outdir, rpt_label, metric, permute)
    #
    return command

def execute(command):
    print command
    os.system(command)

def main(args):
    pair = designParser(args.design)
    command = command_maker(args.jar, args.chip, args.cls, args.exp, args.gmt, args.outdir,\
            args.rpt_label, pair, args.metric, args.permute)
    execute(command)

if __name__=='__main__':
    import os
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--jar', default='/Users/Yoo/YoosScripts/GSEA/gsea2-2.2.4.jar')
    parser.add_argument('-ch', '--chip')
    parser.add_argument('-cl', '--cls')
    parser.add_argument('-e', '--exp')
    parser.add_argument('-g', '--gmt')
    parser.add_argument('-o', '--outdir', default='./')
    parser.add_argument('-rl', '--rpt-label', default='DEG001')
    parser.add_argument('-d', '--design', default='control,case')
    parser.add_argument('-m', '--metric', choices=('Diff_of_Classes', 'Signal2Noise', 'tTest', 'Ratio_of_Classes', 'Euclidean'), default='Signal2Noise')
    parser.add_argument('-p', '--permute', choices=('phenotype', 'geneset'), default='geneset')
    args = parser.parse_args()
    main(args)
