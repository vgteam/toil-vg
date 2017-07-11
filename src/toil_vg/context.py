#!/usr/bin/env python2.7
"""
context.py: Defines a toil_vg context, which contains (and hides) all the config
file values, IOStore dumping parameters, and other things that need to be passed
around but that toil-vg users shouldn't generally need to dinker with.

Instead of dropping the container runner into the command-line options
namespace, we keep them both in here.

"""

import tempfile
import datetime
import pkg_resources

from argparse import Namespace
from toil_vg.vg_config import apply_config_file_args
from toil_vg.vg_common import ContainerRunner, get_container_tool_map
from toil_vg.iostore import IOStore

class Context(object):
    """
    Represents a toil-vg context, necessary to use the library.
    """
    
    def __init__(self, out_store=None, overrides=Namespace()):
        """
        Make a new context, so we can run the library.
        
        Takes an optional Namespace of overrides for default toil-vg
        configuration values.
        
        Overrides can also have a "config" key, in which case that config file
        will be loaded.
        
        """
        
        # Load configuration and apply overrides. If the overrides are from
        # command-line options, we might also get a bunch of tool-specific
        # fields.
        self.config = apply_config_file_args(overrides)
        
        # Make a container runner for running tools
        self.runner = ContainerRunner(container_tool_map=get_container_tool_map(
            self.config))
        
        if out_store is not None:
            # We want to dump files to an IOStore
            # Make it an absolute path while we're getting set up.
            self.out_store_string = IOStore.absolute(out_store)
            
            # Make sure it works by writing some data to it
            f = tempfile.NamedTemporaryFile(delete=True)
            now = datetime.datetime.now()
            f.write('{}\ntoil-vg version {}\nConfiguration:'.format(now,
                pkg_resources.get_distribution('toil-vg').version))
                
            for key, val in self.config.__dict__.items():
                f.write('{}: {}\n'.format(key, val))
            f.flush()
            self.get_out_store().write_output_file(f.name, 'toil-vg-info.txt')
            f.close()
            
        else:
            # We don't want to use an out store
            self.out_store_string = None
            
    def get_out_store(self):
        """
        Return the IOStore to write output files to, or None if they should not
        be written anywhere besides the Toil file store.
        """
        
        if self.out_store_string is not None:
            return IOStore.get(self.out_store_string)
        else:
            return None
            
    def to_options(self, options):
        """
        Turn back into an options-format namespace like toil-vg used to use
        everywhere. Merges with the given real command-line options.
        """
        
        for key, val in self.config.__dict__.items():
            # Copy over everything from the config to the options namespace.
            # Anything in the config that needed to be overridden by command
            # line options already was overridden in the constructor.
            options.__dict__[key] = val

        # Save the other aspects of the context into the options under their old
        # names.
        options.out_store = self.out_store_string
        options.drunner = self.runner
        options.tool = 'library'
        
        return options
        
        
        
        
        
