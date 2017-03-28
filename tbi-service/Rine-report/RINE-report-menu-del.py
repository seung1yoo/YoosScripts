import glob
import os
import shutil

def grep_html_files(dirname, html_files):
    filenames = os.listdir(dirname)
    for filename in filenames:
        full_filename = os.path.join(dirname, filename)
        if os.path.isdir(full_filename):
            grep_html_files(full_filename, html_files)
        else:
            ext = os.path.splitext(full_filename)[-1]
            if ext in ['.html']:
                html_files.append(full_filename)
    return html_files

def modify_html(html_file, menus):
    temp_file = '{0}.tmp'.format(html_file)
    shutil.copy(html_file, temp_file)
    out = open(html_file, 'w')
    for line_ori in open(temp_file):
        line = line_ori.strip()
        if line.startswith('<li class=') and line.endswith('</a></li>'):
            menuname = line.split('</a></li>')[0].split('>')[-1]
            if menuname in menus:
                continue
            #print line_ori
            out.write(line_ori)
        else:
            #print line_ori
            out.write(line_ori)
    out.close()
    os.system('rm {0}'.format(temp_file))


def main(args):
    print(args)
    html_files = []
    html_files = grep_html_files(args.path, html_files)
    for html_file in html_files:
        #print(html_file)
        modify_html(html_file, args.menus)
    print('#DONE')

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', help='RINE report path')
    parser.add_argument('-m', '--menus', nargs='+', choices=('Home','Overview','Samples','Expression','DEG','DEG(Symbol)','GO','Genes','Help'), help='MENU name to delete')
    args = parser.parse_args()
    main(args)
