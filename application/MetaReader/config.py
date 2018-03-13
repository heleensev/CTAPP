import cerberus
import yaml
import os
import sys

predata_schema = "/hpc/local/CentOS7/dhl_ec/software/ctapp/configs/schemas/pre_schema.yaml"
postdata_schema = "/hpc/local/CentOS7/dhl_ec/software/ctapp/configs/schemas/post_schema.yaml"
param_schema = "/hpc/local/CentOS7/dhl_ec/software/ctapp/configs/schemas/param_schema.yaml"
error_log = "/hpc/local/CentOS7/dhl_ec/software/ctapp/logs/preprocess/validator.log"


# class containing info from meta data file as instance attributes
# attributes without setter decorators are checked in ConfigValidator
class Study:

    def __init__(self, id):
        # file properties, manually added in config
        self.studyID = id
        self.path = str
        self.effect_allele = str
        # meta data, added manually in the config
        self.n_studies = int
        self.odds_score = str
        self.phenotype = str
        self.ethnicity = str
        self.build = int
        # file properties of the study, generated in runtime
        self._sep = str
        self._headers = list
        self._head_idx = list
        self._head_to_ref = list
        self._org_headers = list
        self._success = True
        self._disposed = list()
        self._err_mssg = str()
        self._missing = list()

    # update instance attributes directly
    def config_update(self, doc):
        self.__dict__.update(doc)

    @property
    def sep(self):
        return self._sep

    @sep.setter
    def sep(self, s):
        self._sep = s

    @property
    def headers(self):
        return self._headers

    @headers.setter
    def headers(self, heads):
        if not isinstance(heads, list):
            raise Exception('wrong type assignment for "headers')
        self._headers = heads

    @property
    def head_to_ref(self):
        return self._head_to_ref

    @head_to_ref.setter
    def head_to_ref(self, idx_to_ref):
        if not isinstance(idx_to_ref, list):
            raise Exception('wrong type assignment for "head_to_ref')
        self._head_to_ref = idx_to_ref

    @property
    def head_idx(self):
        return self._head_idx

    @head_idx.setter
    def head_idx(self, head_idx):
        if not isinstance(head_idx, list):
            raise Exception('wrong type assignment for "head_idx')
        self._head_idx = head_idx

    @property
    def success(self):
        return self._success

    @success.setter
    def success(self, boolean):
        if not isinstance(boolean, bool):
            raise Exception('wrong type assignment for "success"')
        self._success = boolean

    @property
    def disposed(self):
        return self._disposed

    @disposed.setter
    def disposed(self, dispose):
        if not isinstance(dispose, str):
            raise Exception('wrong type assignment for "disposed"')
        self._disposed.append(dispose)

    @property
    def missing(self):
        return self._missing

    @missing.setter
    def missing(self, missed):
        if not isinstance(missed, list):
            raise Exception('wrong type assignment of "missing')
        self._missing = missed

    @property
    def err_mssg(self):
        return self._err_mssg

    @err_mssg.setter
    def err_mssg(self, mssg):
        if not isinstance(mssg, str):
            raise Exception('wrong type assignment for "err_mssg"')
        self._err_mssg += mssg + '\n'


class Parameters:
    def __init__(self):
        pass

    def dict_update(self, doc):
        self.__dict__.update(doc)


class Resource:
    def __init__(self):
        pass

    def config_update(self, doc):
        self.__dict__.update(doc)


class ConfigValidator:
    def __init__(self, config, preprocess):
        self.terminate = False
        self.preprocess = preprocess
        self.paths = list()
        self.ids = list()
        self.config = self.__read_yaml(config)
        self.__check_form()
        # exit if validation goes south
        if self.terminate:
            print('config file is not valid')
            sys.exit(1)

    def __write_errors(self, errs):
        with open(error_log, 'a') as log:
            log.write(errs)
        self.terminate = True

    def __read_yaml(self, path):
        try:
            with open(path, 'r') as ymlfile:
                ymlconfig = yaml.load(ymlfile)
            return ymlconfig
        except yaml.YAMLError as exc:
            raise Exception("Not a valid YAML: {}".format_map(exc))
        except FileNotFoundError:
            raise Exception("No such file: {}".format(path))

    def __get_schema(self, path):
        schema = self.__read_yaml(path)
        return schema

    def __check_form(self):
        try:
            study_docs = self.config.get('data_info')
            param_docs = self.config.get('parameters')
            for doc in study_docs:
                self.__validate_study(doc)
            for doc in param_docs:
                self.__validate_param(doc)
            self.__check_duplicates('path', self.paths)
            self.__check_duplicates('ID', self.ids)
        except ValueError:
            raise Exception('First level keywords in config not correct')
        except Exception as e:
            raise Exception(e)

    def __validate_param(self, doc):
        schema = self.__get_schema(param_schema)
        v = cerberus.Validator(schema)
        v.validate(doc)
        if v.errors:
            self.__write_errors(yaml.dump(v.errors) + '\n')
        if doc.get('method') == 'subparallel':
            if not doc.get('subchunk'):
                raise Exception('subchunk size should be provided when'
                                'method is "subparallel"')

    def __validate_study(self, doc):
        if self.preprocess:
            data_schema = predata_schema
        else:
            data_schema = postdata_schema
        schema = self.__get_schema(data_schema)
        v = cerberus.Validator(schema)

        ID = doc.get('studyID')
        path = doc.get('meta_data').get('path')
        v.validate(doc)
        if not self.__check_paths(path):
            self.__write_errors('invalid or missing path for: {} \n'.format(ID))
            self.terminate = True
        self.paths.append(path)
        self.ids.append(ID)

        if v.errors:
            self.__write_errors('Errors for study: {}\n'.format(ID))
            self.__write_errors(yaml.dump(v.errors) + '\n')

    def __check_paths(self, p):
        try:
            if os.path.exists(p):
                return True
        except TypeError:
            return False

    def __check_duplicates(self, attr, keys):
        dups = [x for n, x in enumerate(keys) if x in keys[:n]]
        if None in set(dups):
            dups.remove(None)
        if set(dups):
            self.terminate = True
            print('dups: {}'.format(dups))
            self.__write_errors('duplicate {}(s) found: {}'.format(attr, ('\n'.join(dups))))
