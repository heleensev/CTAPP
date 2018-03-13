import sys
import yaml
import gzip
import logging
sys.path.append('/home/sevvy/Scripts/deploy/application')
from application.MetaReader.config import Study, Parameters, Resource

logger = logging.getLogger(__name__)
classes = {'study': Study, 'resource': Resource}


def open_up(fn):
    try:
        with gzip.open(fn) as gzfn:
            gzfn.read(2)
        return gzip.open(fn)
    except OSError:
        info = str(sys.exc_info())
        if 'gzipped' in info:
            return open(fn)
        elif 'FileNotFound' in info:
            raise Exception("No such file: {} in".format(fn, __name__))


def read_meta(path):
    try:
        with open_up(path) as yaml_file:
            meta_yaml = yaml.load(yaml_file)
        return meta_yaml
    except yaml.YAMLError as exc:
        raise Exception("Not a valid YAML: {}".format_map(exc))


def meta_data(path, ob, ix=None):
    objects = list()
    try:
        cls = classes.get(ob)
        config_doc = read_meta(path)
        data_info = config_doc.get('data_info')
        for i, doc in enumerate(data_info):
            if ob == 'study':
                id = doc.get('studyID')
                this_doc = cls(id)
            else:
                this_doc = cls()
            for attrs in ['properties', 'meta_data']:
                subdoc = doc.get(attrs)
                this_doc.config_update(subdoc)
            objects.append(this_doc)
        if isinstance(ix, int):
            objects = objects[ix]
        return objects
    except Exception:
        print(sys.exc_info())


def meta_params(path):
    try:
        config_doc = read_meta(path)
        # the runtime, vmem must be the same for all
        # indices, so only the first doc is fetched
        params = config_doc.get('parameters')[0]
        this_param = Parameters()
        this_param.dict_update(params)
        return this_param
    except Exception:
        print(sys.exc_info())
