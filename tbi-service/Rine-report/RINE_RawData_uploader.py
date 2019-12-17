
import os
import sys

class SAMPLEIN_PARSER:
    def __init__(self, args):
        self.samplein_fn = args.samplein
        #
        self.path = self.select_path()
        self.input_dic = self.select_input()
        #

    def select_path(self):
        for line in open(self.samplein_fn):
            items = line.rstrip('\n').split()
            if items[0] in ['PATH']:
                path = items[1]
            else:
                continue
        return path

    def select_input(self):
        input_dic = dict()
        for line in open(self.samplein_fn):
            items = line.rstrip('\n').split()
            if items[0] in ['INPUT']:
                s_id = int(items[1])
                p_id = int(items[2])
                input_dic.setdefault(s_id,{}).setdefault(p_id,{})
                #
                f_path = items[4].replace('[PATH]',self.path)
                if not os.path.exists(f_path):
                    print 'Not exists : {0}'.format(f_path)
                    sys.exit()
                s_name = items[5]
                input_dic[s_id][p_id].setdefault('path',f_path)
                input_dic[s_id][p_id].setdefault('name',s_name)
                #
            else:
                continue
        return input_dic

class UPLOAD:
    def __init__(self, args, input_dic):
        self.host = args.host
        self.id = args.id
        self.password = args.password
        self.upload_dir = args.upload_dir
        #
        self.link_dir = args.link_dir
        if not os.path.exists(self.link_dir):
            os.makedirs(self.link_dir)
        #
        self.input_dic = input_dic
        #
        self.sym_link_raw()
        self.sym_link_fastqc()

    def sym_link_raw(self):
        for s_id, p_id_dic in sorted(self.input_dic.iteritems()):
            for p_id, info_dic in sorted(p_id_dic.iteritems()):
                #
                target_dir = os.path.join(self.link_dir, info_dic['name'])
                if not os.path.exists(target_dir):
                    os.makedirs(target_dir)
                target_raw = os.path.join(target_dir, '{0}_{1}.fq.gz'.format(info_dic['name'],p_id))
                if not os.path.exists(target_raw):
                    os.symlink(info_dic['path'], target_raw)
                #
    def sym_link_fastqc(self):
        for s_id, p_id_dic in sorted(self.input_dic.iteritems()):
            for p_id, info_dic in sorted(p_id_dic.iteritems()):
                #
                target_dir = os.path.join(self.link_dir, info_dic['name'])
                #
                base_name = info_dic['path'].rstrip('.fq.gz')
                if os.path.exists('{0}_fastqc.zip'.format(base_name)):
                    target = os.path.join(target_dir, '{0}_{1}_fastqc.zip'.format(info_dic['name'],p_id))
                    if not os.path.exists(target): os.symlink('{0}_fastqc.zip'.format(base_name), target)
                if os.path.exists('{0}_fastqc.html'.format(base_name)):
                    target = os.path.join(target_dir, '{0}_{1}_fastqc.html'.format(info_dic['name'],p_id))
                    if not os.path.exists(target): os.symlink('{0}_fastqc.html'.format(base_name), target)
                if os.path.exists('{0}_fastqc'.format(base_name)):
                    target = os.path.join(target_dir, '{0}_{1}_fastqc'.format(info_dic['name'],p_id))
                    if not os.path.exists(target): os.symlink('{0}_fastqc'.format(base_name), target)

    def exe_sshpass(self):
        cmd = "sshpass -p '{0}' scp -P {5} -r {1}/* {2}@{3}:{4}".format(
            self.password, self.link_dir, self.id, self.host, self.upload_dir, args.port)
        print cmd
        os.system(cmd)

def main(args):
    sp = SAMPLEIN_PARSER(args)
    for s_id, p_id_dic in sorted(sp.input_dic.iteritems()):
        for p_id, info_dic in sorted(p_id_dic.iteritems()):
            print s_id, p_id, info_dic['name'], info_dic['path']
    #
    upload = UPLOAD(args, sp.input_dic)
    upload.exe_sshpass()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--samplein', default='./sample.quant.in')
    parser.add_argument('--link-dir', default='./Upload')
    parser.add_argument('--host', default='59.18.159.30')
    parser.add_argument('--port', default='2022')
    parser.add_argument('--id', default='CDC_Yjh')
    parser.add_argument('--password', default='CDC_dlvmflxm2YJH')
    parser.add_argument('--upload-dir', default='/BiO/Live/CDC_Yjh/TBD180495_20180827/0.RawData')
    args = parser.parse_args()
    main(args)
