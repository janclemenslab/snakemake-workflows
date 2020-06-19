import sys
try:
    import PySide2  # this will force pyqtgraph to use PySide instead of PyQt4/5
except ImportError:
    pass

import tkinter.filedialog as filedialog


from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import formbuilder
from pprint import pprint
from collections import defaultdict
import logging
import yaml
import os

from videoreader import VideoReader
import os

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
    return rules

# movie_name = f"dat/localhost-20200110_094433/localhost-20200110_094433.mp4"
# movie_name = 'Z:/#Common/chainingmic/dat/localhost-20200617_185618/localhost-20200617_185618.mp4'#os.path.abspath(movie_name)


# do the damage
if len(sys.argv) > 1:
    movieName = sys.argv[1]
else:
    movieName = filedialog.askopenfilename()
movie_name = movieName

print(movie_name)
basename = os.path.basename(os.path.splitext(movie_name)[0])
print(basename)
rootname = movie_name.split('dat')[0][:-1]
print(rootname)
dirname = os.path.basename(rootname)
print(dirname)
localdir = os.path.dirname(os.path.realpath(__file__))

app = QtGui.QApplication([])

yaml_file = os.path.join(localdir, "analysis profiles.yaml")
with open(yaml_file, "r") as form_yaml:
    items_to_create = yaml.load(form_yaml, Loader=Loader)

# autopopulate list with profiles from folder:
profile_dir = os.path.join(rootname, 'workflow', 'analysis_profiles')

# look in `workflow/analysis profiles` in the folder of the current experiment
# also include defaults from `snakemake-workflows/analysis profiles` (at least new.yaml)
profiles = os.listdir(profile_dir)

for profile in profiles:
    profile_name = profile#os.path.splitext(profile)[0]
    profile_path = os.path.join(profile_dir, profile)
    # append to pull down options
    items_to_create['main'][0]['options'] += "," + profile_name
    # create sub-form
    with open(profile_path, "r") as form_yaml:
        sub_form = yaml.load(form_yaml, Loader=Loader)
    items_to_create['main'][0][profile_name] = sub_form

# find a way to specify defaults for each folder - maybe instead of "None",
# set it to a specified "analysis profiles/default.yaml" if that file exists
if 'default.yaml' in profiles:
    items_to_create['main'][0]['default'] = 'default.yaml'

form = formbuilder.DictFormWidget(form_dict=items_to_create)

dialog = formbuilder.FormBuilderModalDialog(form, None)
data = dialog.get_results()

if data is None:
    logging.info('will not add custom profile - this means it will use (DESCRIBE DEFAULTS HERE!')
else:
    logging.debug(data)
    pprint(data)
    rules = dot_keys_to_nested(data)
    pprint(dict(rules))
    # TODO: clean up paths to be relative
    logging.info(dict(rules))
    filename = os.path.splitext(movie_name)[0] + '_analysis.yaml'
    logging.info(f"saving to {filename}")
    with open(filename, 'w') as f:
        yaml.dump(dict(rules), f)
