import yaml
import os
import shutil
from collections import defaultdict
from pprint import pprint


import yaml
from typing import Any, IO, Dict

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


def read_annotations(wildcards, input):
    """
    input (List[str]) - format `stuff/dirname/stuff`
    """
    dirname = os.path.basename(os.path.dirname(input[0]))
    param_file = f'dat/{dirname}/{dirname}_analysis.yaml'
    if not os.path.exists(param_file):
        # get defailts from default schema
        schema_file = 'workflow/analysis_profiles/default.yaml'
        with open(schema_file, 'r') as fid:
            schema = yaml.load(fid, Loader=Loader)
        params = params_from_schema_defaults(schema)
        params = dot_keys_to_nested(params)
    else:
        with open(param_file, 'r') as fid:
            params = yaml.load(fid, Loader=yaml.FullLoader)
    return params

def get_all_targets_for_expt(wildcards):
    targets = defaultdict(lambda: [])
    for wildcard in wildcards:
        inp = f"dat/{wildcard}/{wildcard}.ext"
        analysis = read_annotations(None, [inp])
        for key in analysis.keys():
            suffix = analysis[key]['suffix'] if 'suffix' in analysis[key] else key
            extension = analysis[key]['extension'] if 'extension' in analysis[key] else 'h5'
            targets[key].append(f"res/{wildcard}/{wildcard}_{suffix}.{extension}")
    return targets


## ALL OF THE ABOVE SHOULD GO TO snake-workflows/utils.py

# list all video files
directories, files, = glob_wildcards('dat/{directory}/{file}.mp4')

# remove hidden
directories = [directory for directory, file in zip(directories, files) if not file.startswith('.')]
files = [file for directory, file in zip(directories, files) if not file.startswith('.')]

# check analysis.yaml for all files and decide which of these
# to add to the output file list in the next (based on existence of a block with the rule name)
targets = defaultdict(lambda: [])
for dire in directories:
    analysis = read_annotations(None, [f'dat/{dire}/{dire}',])
    for key in analysis.keys():
        extension = analysis[key]['extension'] if 'extension' in analysis[key] else 'h5'
        suffix = analysis[key]['suffix'] if 'suffix' in analysis[key] else key
        targets[key].append(f"res/{dire}/{dire}_{suffix}.{extension}")
        # add all to the clean rule
    targets['clean'].append(f"res/{dire}/DONE")
    targets['move'].append(f"dat.processed/{dire}")

pprint(dict(targets))