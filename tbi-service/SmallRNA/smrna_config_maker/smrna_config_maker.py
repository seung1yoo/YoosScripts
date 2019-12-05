

class SMRNA_CONFIG_MAKER:
    conf_temp_header = '''<?xml version='1.0' encoding='UTF-8' ?>\n<CONFIG>\n'''
    conf_temp_parameters = '''\
    <PARAMETERS>
        <PROJECT>[PROJECT]</PROJECT>
        <PROJECT_ID>[PROJECT]</PROJECT_ID>
        <PROJECT_DIR>[PATH]/Result</PROJECT_DIR>
    </PARAMETERS>
'''
    conf_temp_data_open = '''\
    <DATA>
'''
    conf_temp_data = '''\
        <FASTQ ID='[DATA_ID]' SM='[DATA_SM]' PL='ILLUMINA' SCALE='33' KIT='nextflex'>
            <FORWARD>[RAW_PATH]</FORWARD>
        </FASTQ>
'''
    conf_temp_data_close = '''\
    </DATA>
'''
    conf_temp_deg = '''\
    <DEG name='[DEG_ID]'>
        [DEG_ID]:SAMPLES=[DEG_SAMPLES]
        [DEG_ID]:GROUPS=[DEG_GROUPS]
    </DEG>
'''
    conf_temp_footer = '''</CONFIG>\n'''

    def __init__(self, args):
        self.project = args.project
        self.path = args.path
        self.maker_conf = args.maker_conf
        self.conf_dic = self.read_maker_conf()
        #
        self.out_fh = open(args.out, 'w')
        self.out_fh.write(self.conf_temp_header)
        conf_parameters = self.conf_temp_parameters
        conf_parameters = conf_parameters.replace('[PROJECT]', self.project)
        conf_parameters = conf_parameters.replace('[PATH]', self.path)
        self.out_fh.write(conf_parameters)

        self.out_fh.write(self.conf_temp_data_open)
        for data_sm, info_s in sorted(self.conf_dic['INPUT'].items(), key=lambda (k,v): (v,k)):
            conf_data = self.conf_temp_data
            conf_data = conf_data.replace('[DATA_ID]', info_s[0])
            conf_data = conf_data.replace('[DATA_SM]', data_sm)
            conf_data = conf_data.replace('[RAW_PATH]', info_s[1])
            self.out_fh.write(conf_data)
        self.out_fh.write(self.conf_temp_data_close)

        for deg_id, info_dic in sorted(self.conf_dic['DEG'].items()):
            conf_deg = self.conf_temp_deg
            conf_deg = conf_deg.replace('[DEG_ID]', deg_id)
            conf_deg = conf_deg.replace('[DEG_SAMPLES]', info_dic['samples'])
            conf_deg = conf_deg.replace('[DEG_GROUPS]', info_dic['groups'])
            self.out_fh.write(conf_deg)

        self.out_fh.write(self.conf_temp_footer)
        self.out_fh.close()

    def read_maker_conf(self):
        data_num = 0
        deg_num  = 0
        conf_dic = dict()
        for line in open(self.maker_conf):
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split('\t')
            conf_tag = items[0]
            if conf_tag in ['INPUT']:
                data_num += 1
                data_id = 'S{0:0>4}'.format(data_num)
                data_sm = items[1]
                data_path = items[2]
                conf_dic.setdefault(conf_tag, {}).setdefault(data_sm, (data_id, data_path))
            elif items[0] in ['DEG']:
                deg_num += 1
                deg_id = 'DEG_{0:0>3}'.format(deg_num)
                ctrl_samples = [conf_dic['INPUT'][data_sm][0] for data_sm in items[1].split(':')]
                case_samples = [conf_dic['INPUT'][data_sm][0] for data_sm in items[2].split(':')]
                ctrl_groups  = ['Control' for data_sm in items[1].split(':')]
                case_groups  = ['Case' for data_sm in items[2].split(':')]
                conf_dic.setdefault(conf_tag, {}).setdefault(deg_id, {}).setdefault(
                        'samples', '{0}:{1}'.format(':'.join(ctrl_samples), ':'.join(case_samples)))
                conf_dic.setdefault(conf_tag, {}).setdefault(deg_id, {}).setdefault(
                        'groups', '{0}:{1}'.format(':'.join(ctrl_groups), ':'.join(case_groups)))
            else:
                print('error')
                import sys
                sys.exit()
        return conf_dic


def main(args):
    maker = SMRNA_CONFIG_MAKER(args)
    print(maker.conf_dic)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--project',
            default='TBD000000')
    parser.add_argument('--path',
            default='/BiO/BioPeople/siyoo/CNU-Cow-smallRNA-20190424')
    parser.add_argument('--maker-conf',
            default='/BiO/BioPeople/siyoo/CNU-Cow-smallRNA-20190424/Raw/smrna_config_maker.conf')
    parser.add_argument('--out',
            default='config.xml')
    args = parser.parse_args()
    main(args)
