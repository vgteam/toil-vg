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
    add_common_vg_parse_args, add_container_tool_parse_args
from toil_vg.vg_mapeval import run_map_eval_summarize, run_map_eval_table, run_map_eval_plot
from toil_vg.vg_calleval import run_calleval_plots
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
    
    # Mapeval plotting options
    parser.add_argument('--position-stats', type=make_url, default=None,
                        help='position.results.tsv file from a mapeval run')
    parser.add_argument('--tables-only', action='store_true',
                        help='make only summary tables and not plots')
    parser.add_argument('--plots-only', action='store_true',
                        help='make only plots and not summary tables')
    
    # Calleval plotting options
    parser.add_argument('--roc-base',
                        help='ROC file base URL, under which <condition>_vcfeval_output/<type>_roc.tsv.gz files exist')
    parser.add_argument('--names', nargs='+',
                        help='condition names for ROC files to load')
    parser.add_argument('--clipping', nargs='+', default=['clipped'],
                        help='clipping modes to load')
    
    # General options
    parser.add_argument('--plot-sets', nargs='+', default=[],
                        help='comma-separated lists of condition-tagged GAM names (primary-mp-pe, etc.) with colon-separated title prefixes')
    
                        
    # We also need to have these options to make lower-level toil-vg code happy
    # with the options namespace we hand it.
    
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)
    
def validate_options(options):
    """
    Throw an error if an invalid combination of options has been selected.
    """

    
    require(options.position_stats or options.roc_base, 'either --position-stats or --roc-base is required')
    require(not (options.position_stats and options.roc_base), 'cannot operate on --position-stats and --roc-base in the same run')
    require(options.roc_base is None or options.names is not None, '--roc-base requires --names')
    require(options.names is None or options.roc_base is not None, '--names requires --roc-base')
    
    
def run_plot(job, context, options, position_stats_file_id=None, eval_results_dict=None, plot_sets=[None]):
    """
    Main Toil job.
    
    If position_stats_file_id is given, make mapping plots and/or tables based
    on the given position_stats_file_id stats file, restricting to the given
    set of sets of conditions to plot.
    
    If instead eval_results_dict is given, as a dict from condition name, then
    'clipped'/'unclipped', then 'snp'/'non_snp'/'weighted' to ROC data file ID,
    then plot the variant calling ROC curves.
    
    """
    
    if position_stats_file_id is not None:
        # Do position stats
    
        # Do plots and tables, unless we want to do just tables
        job_fn = run_map_eval_summarize
        if options.tables_only:
            job_fn = run_map_eval_table
        if options.plots_only:
            job_fn = run_map_eval_plot
        
        plot_job = job.addChildJobFn(job_fn, context, position_stats_file_id, plot_sets,
                                     cores=context.config.misc_cores, memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
                                     
    elif eval_results_dict is not None:
        # Do variant calling plotting
        
        # Get the names back
        names = eval_results_dict.keys()
        
        plot_job = job.addChildJobFn(run_calleval_plots, context, names, eval_results_dict, plot_sets=plot_sets)
    
    else:
        raise RuntimeError('No position stats or vcfeval results available!')

def parse_plot_set(plot_set_string):
    """
    
    Given one of the string arguments to the --plot-sets option, parse out a
    data structure representing which conditions ought to be compared against
    each other, and what those comparison plots/tables should be called.

    The syntax of a plot set is [title:]condition[,condition[,condition...]].
    
    The first condition is the comparison baseline, when applicable.
    
    Returns a tuple of a plot set title, or None if unspecified, and a list of
    condition names.
    
    """
    
    colon_pos = plot_set_string.find(':')
    
    if colon_pos != -1:
        # Pull out the title before the colon
        title = plot_set_string[0:colon_pos]
        # And the rest of the specifier after it
        plot_set_string = plot_set_string[colon_pos + 1:]
    else:
        # No title given
        title = None
        
    # Return the title and condition list tuple
    return (title, plot_set_string.split(','))
        
    
def parse_plot_sets(plot_sets_list):
    """
    
    Given a list of plot set strings, parses each with parse_plot_set.
    
    Returns a list of tuples. Each tuple is a plot set title, or None if no
    title is to be applied, and a list of condition names, or None if all
    conditions are to be included.
    
    If no plot sets are specified in the list, returns a single plot set for
    all conditions.
    
    """
    
    plot_sets = [parse_plot_set(spec) for spec in plot_sets_list]
    if len(plot_sets) == 0:
        # We want to plot everything together
        # We use the special None value to request that.
        plot_sets = [('All', None)]
    
    return plot_sets
    
def title_to_filename(kind, i, title, extension):
    """
    Given the kind of thing you want to save ('table', 'plot-qq', etc.), the
    number of the thign out of all things of that type, the human-readable
    title ('All Conditions vs. Whatever'), and an extansion, come up with a
    safe filename to save the plot under.
    
    The title may be None, in which case it is ommitted.
    
    The extension may be None, in which case it is omitted.
    """
    
    if title is not None:
        # Filter down to good filename characters as in https://stackoverflow.com/a/7406369
        safe_title = ''.join((c for c in title if c.isalnum()))
    else:
        safe_title = None
        
    # The name always includes the kind of thing
    part_list = [kind]
    
    if i != 0:
        # Include number padded to 2 digits only if nonzero. If not included,
        # the next dash will sort before the other numbers. TODO: Add the 00 in
        # here and break backward compatibility with everything looking for the
        # output files.
        part_list.append('-{:02d}'.format(i))
        
    if title is not None:
        # Filter down to good filename characters as in https://stackoverflow.com/a/7406369
        safe_title = ''.join((c for c in title if c.isalnum()))
        part_list.append('-{}'.format(safe_title))
        
    if extension is not None:
        part_list.append('.{}'.format(extension))
        
    return ''.join(part_list)

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
    
    plan.position_stats_file_id = None
    if options.position_stats:
        # We will use a position stats file with all the reads to re-plot and
        # their correctness statuses and MAPQs
        plan.position_stats_file_id = toil.importFile(options.position_stats)
    
    plan.eval_results_dict = None
    if options.roc_base:
        # We will load up a bunch of ROCs
        
        plan.eval_results_dict = collections.defaultdict(lambda: collections.defaultdict(dict))
        
        for condition in options.names:
            # For each condition that should have a ROC
            
            for clipping in options.clipping:
                # For each clipping mode
                
                # What should the condition be with the clip tag?
                # If we are running just one mode it is just the condition.
                # But if we are running both it is the condition and the condition-unclipped
                clip_tag_condition = condition
                if clipping == 'unclipped' and 'clipped' in options.clipping:
                    clip_tag_condition += '-unclipped'
                    
                for roc_name in ['snp', 'non_snp', 'weighted']:
                
                    # What URL do we want
                    url = make_url(options.roc_base + '/' + clip_tag_condition + '_vcfeval_output/' + roc_name + '_roc.tsv.gz')
                    
                    # Load it up in the appropriate spot
                    plan.eval_results_dict[condition][clipping][roc_name] = toil.importFile(url)
   
    # Also process options that need parsing
    plan.plot_sets = parse_plot_sets(options.plot_sets)
   
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
                                     position_stats_file_id=plan.position_stats_file_id,
                                     eval_results_dict=plan.eval_results_dict,
                                     plot_sets=plan.plot_sets)
                
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
    
