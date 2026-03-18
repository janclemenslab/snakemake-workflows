import yaml
import os
import shutil
from collections import defaultdict
import yaml
from typing import Any, IO, Dict
import logging


class Loader(yaml.SafeLoader):
    """YAML Loader with `!include` constructor."""

    def __init__(self, stream: IO) -> None:
        """Initialise Loader."""

        try:
            self._root = os.path.split(stream.name)[0]
        except AttributeError:
            self._root = os.path.curdir

        super().__init__(stream)


def construct_include(loader: Loader, node: yaml.Node) -> Any:
    """Include file referenced at node."""

    filename = os.path.abspath(os.path.join(loader._root, loader.construct_scalar(node)))
    extension = os.path.splitext(filename)[1].lstrip('.')

    with open(filename, 'r') as f:
        if extension in ('yaml', 'yml'):
            return yaml.load(f, Loader)
        else:
            return ''.join(f.readlines())


yaml.add_constructor('!include', construct_include, Loader)


def dot_keys_to_nested(data: Dict) -> Dict:
    """old['aaaa.bbbb'] -> d['aaaa']['bbbb']

    Args:
        data (Dict): [description]

    Returns:
        Dict: [description]
    """
    rules = defaultdict(lambda: dict())
    for key, val in data.items():
        if '.' in key:
            key, _, param = key.partition('.')
            rules[key][param] = val
    return dict(rules)  # make sure we return a normal dict - defaultdict is not pickle-able


def params_from_schema_defaults(schema: Dict) -> Dict:
    default_params = dict()
    for entry in schema:
        try:
            default_params[entry['name']] = entry['default']
        except KeyError:
            pass
    return default_params


def read_annotations_boxes(wildcards, input):
    """
    input (List[str]) - format `stuff/dirname/stuff`
    """
    print("wildcards:", wildcards, "input:", input)
    dirname = os.path.basename(os.path.dirname(input[0]))[:-6]
    param_file = f'dat/{dirname}/{dirname}_analysis.yaml'
    if not os.path.exists(param_file):
        # print(param_file, ' not found')
        # get defailts from default schema
        schema_file = 'workflow/analysis_profiles/default_sleap_boxes.yaml'
        with open(schema_file, 'r') as fid:
            schema = yaml.load(fid, Loader=Loader)
        params = params_from_schema_defaults(schema)
        params = dot_keys_to_nested(params)
    else:
        with open(param_file, 'r') as fid:
            params = yaml.load(fid, Loader=yaml.FullLoader)
    return params

def read_annotations(wildcards, input):
    """
    input (List[str]) - format `stuff/dirname/stuff`
    """
    # print("wildcards:", wildcards, "input:", input)
    # print("os.path.basename(os.path.dirname(input[0]))", os.path.basename(os.path.dirname(input[0])))
    # print("os.path.splitext(os.path.split(input[0])[1])[0]", os.path.splitext(os.path.split(input[0])[1])[0])

    # try simple name pattern
    dirname = os.path.basename(os.path.dirname(input[0]))
    param_file = f'dat/{dirname}/{dirname}_analysis.yaml'

    # try alternative name pattern
    if not os.path.exists(param_file):
        dirname = os.path.splitext(os.path.split(input[0])[1])[0]  # this should work for all whatever/{dirname}/whatever
        param_file = f'dat/{dirname}/{dirname}_analysis.yaml'

    if not os.path.exists(param_file):
        # print(param_file, 'NOT found')
        # get defailts from default schema
        schema_file = 'workflow/analysis_profiles/default.yaml'
        with open(schema_file, 'r') as fid:
            schema = yaml.load(fid, Loader=Loader)
        params = params_from_schema_defaults(schema)
        params = dot_keys_to_nested(params)
    else:
        # print(param_file, 'found')
        with open(param_file, 'r') as fid:
            params = yaml.load(fid, Loader=yaml.FullLoader)
        # print(params)

    # detect and pre-process new-style analysis file
    is_new_style =  'Jobs' in params
    if is_new_style:
        params = params['Jobs']

    return params

def get_all_targets_for_expt(wildcards):
    targets = defaultdict(lambda: [])
    directories = [wildcards.dir] if hasattr(wildcards, 'dir') else [wildcards]

    # Keep this in sync with enabled jobs in project workflows.
    enabled_jobs = {
        'tracks': {'suffix': 'tracks', 'extension': 'h5'},
        'sleap': {'suffix': 'sleap', 'extension': 'h5'},
        'song': {'suffix': 'annotations', 'extension': 'csv'},
    }

    for wildcard in directories:
        inp = f"dat/{wildcard}/{wildcard}.ext"
        analysis = read_annotations(None, [inp])
        if isinstance(analysis, list):
            analysis = analysis[0] if analysis else {}
        if not isinstance(analysis, dict):
            continue

        for key, defaults in enabled_jobs.items():
            job = analysis.get(key)
            if not isinstance(job, dict):
                continue

            suffix = job.get('suffix', defaults['suffix'])
            extension = job.get('extension', defaults['extension'])
            if extension:
                targets[key].append(f"res/{wildcard}/{wildcard}_{suffix}.{extension}")
            else:
                targets[key].append(f"res/{wildcard}/{wildcard}_{suffix}")
    return targets


## ALL OF THE ABOVE SHOULD GO TO snake-workflows/utils.py

# list all video files
directories, files, = glob_wildcards('dat/{dir}/{file}.mp4')

# remove hidden
directories = [dir for dir, file in zip(directories, files) if not file.startswith('.')]
# files = [file for dir, file in zip(directories, files) if not file.startswith('.')]
# import rich
# check analysis.yaml for all files and decide which of these
# to add to the output file list in the next (based on existence of a block with the rule name)
targets = defaultdict(lambda: [])

for dir in directories:
    logging.debug(f'dat/{dir}/{dir}')
    analysis = read_annotations(None, [f'dat/{dir}/{dir}',])
    logging.debug(analysis)
    if isinstance(analysis, list):
        analysis = analysis[0]
    for key in analysis.keys():
        extension = analysis[key]['extension'] if 'extension' in analysis[key] else 'h5'
        suffix = analysis[key]['suffix'] if 'suffix' in analysis[key] else key
        if len(extension):
            targets[key].append(f"res/{dir}/{dir}_{suffix}.{extension}")
        else:
            targets[key].append(f"res/{dir}/{dir}_{suffix}")
        # add all to the clean rule
    targets['clean'].append(f"res/{dir}/DONE")
    targets['move'].append(f"dat.processed/{dir}")
logging.debug(dir)
