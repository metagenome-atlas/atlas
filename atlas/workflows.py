import os
from atlas.utils import validate_assembly_config


def assemble(config):
    validate_assembly_config(config)
    os.path.join(os.path.dirname(os.path.abspath(__file__))
