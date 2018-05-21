#!/usr/bin/env python2.7
"""
vg_plot.py: Make plots from toil-vg experiments

"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string, math
import urlparse
import getpass
import pdb
import gzip
import logging
import copy
from collections import Counter

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import require, make_url, remove_ext,\
    add_common_vg_parse_args, add_container_tool_parse_args, get_vg_script
from toil_vg.vg_mapeval import run_map_eval_summarize, run_map_eval_table
from toil_vg.context import Context, run_write_info_to_outstore

logger = logging.getLogger(__name__)

def plot_subparser(parser):
    """
    Create a subparser for plot.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # Add the out_store
    # TODO: do this at a higher level?
    # Or roll into Context?
    parser.add_argument('out_store',
                        help='output store.  All output written here. Path specified using same syntax as toil jobStore')
    
    # Add plot stuff
    add_plot_options(parser)
    
    # Add common docker options
    add_container_tool_parse_args(parser)
    
def add_plot_options(parser):
    """
    Add the mapeval options to the given argparse parser.
    """
    
    # General options
    parser.add_argument('--position-stats', type=make_url, default=None,
                        help='position.results.tsv file from a mapeval run')
    parser.add_argument('--plot-sets', nargs='+', default=[],
                        help='comma-separated lists of condition-tagged GAM names (primary-mp-pe, etc.) to plot together')
    parser.add_argument('--tables-only', action='store_true',
                        help='make only summary tables and not plots')
                        
    # We also need to have these options to make lower-level toil-vg code happy
    # with the options namespace we hand it.
    
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)
    
def validate_options(options):
    """
    Throw an error if an invalid combination of options has been selected.
    """

    require(options.position_stats, 'a --position-stats file is required')
    
def run_plot(job, context, options, position_stats_file_id, plot_sets):
    """
    Main Toil job, and main entrypoint for use as a library.
    
    Plot plots based on the given position_stats_file_id stats file,
    restricting to the given set of sets of conditions to plot.
    
    """
    
    # Do plots and tables, unless we want to do just tables
    job_fn = run_map_eval_summarize if not options.tables_only else run_map_eval_table
    
    plot_job = job.addChildJobFn(job_fn, context, position_stats_file_id, plot_sets,
                                 cores=context.config.misc_cores, memory=context.config.misc_mem,
                                 disk=context.config.misc_disk)

def make_plot_plan(toil, options):
    """
    Import all the necessary files form options into Toil.
    
    Keep the IDs under names in an argparse namespace that functions as a "plan"
    for the workflow.
    
    """
    
    # Make a plan
    plan = argparse.Namespace()
    
    start_time = timeit.default_timer()
            
    # Upload local files to the remote IO Store
    
    if options.position_stats:
        # We will use a position stats file with all the reads to re-plot and
        # their correctness statuses and MAPQs
        plan.position_stats_file_id = toil.importFile(options.position_stats)
    else:
        plan.position_stats_file_id = None
   
    # Also process options that need parsing
    # TODO: Make mapeval include the plot sets in its plan as well?
    # TODO: Pass the plan along instead of getting theings from the options or unpacking it?
   
    plan.plot_sets = [spec.split(',') for spec in options.plot_sets]
    if len(plan.plot_sets) == 0:
        # We want to plot everything together
        # We use the special None value to request that.
        plan.plot_sets = [None]

    end_time = timeit.default_timer()
    logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))
    
    return plan

def plot_main(context, options):
    """
    Run the mapeval workflow.
    """

    # Make sure the options are good
    validate_options(options)
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    t = copy.deepcopy(context)
    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            # Import everything
            plan = make_plot_plan(toil, options)

            # Make a job to run the mapeval workflow, using all these various imported files.
            main_job = Job.wrapJobFn(run_plot,
                                     context, 
                                     options, 
                                     plan.position_stats_file_id,
                                     plan.plot_sets)
                
            # Output files all live in the out_store, but if we wanted to we could export them also/instead.

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv)
            init_job.addFollowOn(main_job)

            # Run the root job
            toil.start(init_job)
        else:
            toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
