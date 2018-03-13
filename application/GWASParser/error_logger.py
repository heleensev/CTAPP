import pandas as pd
from os import path

config_dir = '/home/dhl_ec/hseverin/deploy/application/configs/preprocess'
logdir = '/home/dhl_ec/hseverin/deploy/logs/preprocess'


def init_usr_check(studies):

    failed = [study for study in studies if not study.success]
    succeeded = [study for study in studies if study.success]

    for studies, log_write in zip([succeeded, failed], (succeeded_logs, failed_log)):
        dfs = list()
        for study in studies:
            dfs.append(get_heads(study))
        if dfs:
            log_write(dfs, studies)


def get_heads(study):
    head = pd.read_csv(filepath_or_buffer=study.path, sep=study.sep, nrows=10,
                       names=study.headers, usecols=study.head_idx, header=0)
    return head


def write_heads(df):
    out_file = path.join(logdir, 'preparser.log')
    with open(out_file, 'a') as out_file:
        df.to_csv(path_or_buf=out_file, sep='\t', mode='a')
        out_file.write('\n\n')


def writer(to_write):
    logpath = path.join(logdir, 'preparser.log')
    with open(logpath, 'a') as out:
        out.write(to_write)


def succeeded_logs(dfs, studies):
    writer('The classifying of columns was successful for the following studies (hurah)\n')
    for st, df in zip(studies, dfs):
        print(st.studyID)
        writer('Study: {}\n'.format(st.studyID))
        write_heads(df)
    writer('\n_________________________________________________\n')


def failed_log(dfs, studies):
    writer('The classifying of columns was NOT successful for the following studies:\n')
    for st, df in zip(studies, dfs):
        writer('Study {}\n'.format(st.studyID))
        if st.err_mssg:
            writer('Error messages for study:\n{}'.format(st.err_mssg))
        write_heads(df)
    writer('Please update the config file (in {})\n'
           'Make sure to processor.py again like: processor.py validate --config *config*\n'.format(config_dir))
    writer('\n_________________________________________________\n')


