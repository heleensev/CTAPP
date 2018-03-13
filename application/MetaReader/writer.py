import yaml
import logging
import traceback
import copy
import timebuddy

logger = logging.getLogger(__name__)
date = ''.join(timebuddy.date().split('-'))
updated_config = '/hpc/local/CentOS7/dhl_ec/software/ctapp/configs/processed/update{}.yaml'.format(date)


class ConfigError(Exception):
    # always exit program if error occurs
    def __init__(self, message):
        super(ConfigError, self).__init__(message)


def update_meta(meta_doc, study_list):
    meta_doc['study_info'] = study_list
    return meta_doc


def read_yaml(path):
    try:
        with open(path, 'r') as ymlfile:
            doc = yaml.load(ymlfile)
        return doc
    except yaml.YAMLError as exc:
        raise Exception("Not a valid YAML: {}".format_map(exc))
    except FileNotFoundError:
        raise Exception("No such file: {}".format(path))


def write_config(config):
    try:
        with open(updated_config, 'w') as jfile:
            doc = yaml.dump(config, jfile, default_flow_style=False)
        return doc
    except Exception as e:
        logger.error(traceback.format_exc())
        raise ConfigError(e)


def build_config(studies, param):

    study_temp = {
        'properties': {}
        , 'meta_data': {}
    }
    config_update = dict(data_info=[], parameters=[])

    param_info = config_update['parameters']
    study_info = config_update['data_info']

    param_doc = dict()
    param_doc['chunksize'] = param.chunksize
    param_doc['subchunk'] = param.subchunk
    param_doc['operation'] = param.operation
    param_doc['method'] = param.method
    param_doc['mongo'] = param.mongo
    param_doc['runtime'] = param.runtime
    param_doc['vmem'] = param.vmem

    param_info.append(param_doc)

    for study in studies:
        study_doc = copy.deepcopy(study_temp)
        study_doc['studyID'] = study.studyID

        props = study_doc['properties']
        props['headers'] = study.headers
        props['header_indices'] = study.head_idx
        props['separator'] = study.sep
        props['skip'] = study.skip

        info = study_doc['meta_data']
        info['path'] = study.path
        info['phenotype'] = study.phenotype
        info['effect_type'] = study.effect_type
        # non-essential attributes may be null
        if hasattr(study, 'n_studies'):
            info['n_studies'] = study.n_studies
        else:
            info['n_studies'] = None
        if hasattr(study, 'ethnicity'):
            info['ethnicity'] = study.ethnicity
        else:
            info['ethnicity'] = 'EUR'
        if hasattr(study, 'build'):
            info['build'] = study.build
        else:
            info['build'] = None
        study_info.append(study_doc)

    write_config(config_update)
