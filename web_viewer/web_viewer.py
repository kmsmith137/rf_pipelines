#!/usr/bin/env python

from os import environ, listdir
from os.path import isfile, isdir, exists, join
from json import loads
from math import ceil
from flask import Flask, url_for, send_from_directory
app = Flask(__name__)


"""
This is web viewer for the L1 pipeline. Currently, it retrieves files from the /data2/web_viewer directory, 
but it can be easily modified to get files from multiple machines. Running the web viewer will display
a list of users, each with a list of pipeline runs in their directories.

A persistent web viewer is running from the frb1 web_viewer account, and is up at frb1.physics.mcgill.ca:5000/!
    Show tiles - displays all outputted plots (default: zoom 0*, index1 0, index2 4)
    Show triggers - displays all triggers at a specified zoom - must be modified from url (defult: 0)
    Show last transform - same as show triggers, but for the last plotter transform (note this will display
    trigger plots if the bonsai plotter produced multiple groups of plots (with the multi-tree plotter)

DEPENDENCIES
Flask (pip install Flask)

RUNNING
In your web_viewer directory, 
    ./run-web-viewer.sh /data2/web_viewer 5000

*Slightly annoying note: the zoom levels in the url are opposite those in the filenames, hence the 
 'reverse()' in the Parser class. This is only really relevant if a user is modifying the zoom level
 by hand from the url. 
"""


class Parser():
    """
    This gets the list of file names at different zoom levels produced by the plotter transforms, and reads 
    all information relevant to displaying the required images (min/max index and zoom values and timestamps)
    for an individual pipeline run. This is repeatedly called by the Crawler() class as it discovers new 
    users and pipeline run directories. 

    Members
    -------

      self.fnames[itransform][izoom][it] -> string
      self.ftimes[itransform][izoom][it] -> float
      self.min_zoom -> int
      self.max_zoom -> int
      self.min_index -> int
      self.max_index[itransform][izoom] -> int

    The 'itransform' index runs over all plot_groups (e.g. one or more waterfall plots, followed by one
    or more bonsai trigger plots).  

    The 'izoom' index runs over zoom levels, ordered from most-zoomed-out to most-zoomed-in.  (This ordering
    is opposite to the file naming convention in the rf_pipelines.plotter_transform.)
    """

    def __init__(self, path):
        # First, read everything in /data2/web_viewer, and extract file names and timestamps
        # Initializes self.fnames, self.ftimes.
        self._get_files(path) 

        # Check whether there was anything of value in the run
        if len(self.fnames) == 0:
            self.fnames = self.ftimes = None

        # Calculate some useful values for the viewer, given a successful run
        if self.fnames is not None and self._check_zoom():
            self.min_zoom, self.min_index = 0, 0
            self.max_zoom = len(self.fnames[0])
            self.max_index = [[len(zoom) for zoom in transform] for transform in self.fnames]

        # In case a directory does not contain plots or transforms contain a different number of zoom levels
        else:
            self.min_zoom, self.min_index = None, None
            self.max_zoom = None
            self.max_index = None
        
    def _get_files(self, path):
        """
        Initializes self.fnames (plot filenames) and self.ftimes (plot start times), based on the .json file 
        produced from pipeline runs.  The output is in the following form:
        
        [[[z0tf0f0, z0tf0f1, ...], [z1tf0f0, z1tf0f1, ...], ..., [...]],
         [[z0tf1f0, z0tf1f0, ...], [z1tf1f0, z1tf1f1, ...], ..., [...]],
         [...]]
        """

        json_file = open(path + '/rf_pipeline_0.json').read()
        json_data = loads(json_file)

        self.t0 = json_data['t0']
        self.s_per_sample = (json_data['t1'] - json_data['t0']) / json_data['nsamples']  # number of seconds per sample

        self.fnames = []  # stores file names (so the viewer can request them)
        self.ftimes = []  # stores timestamps (to be displayed for chime_stream_from_times())
        
        for transform in json_data['transforms']:
            if ('plotter_transform' in transform['name']) and ('plots' in transform):
                self._add_plotter_transform(transform)
            elif ('bonsai_dedisperser' in transform['name']) and ('plots' in transform):
                self._add_bonsai_dedisperser(transform)


    def _add_plotter_transform(self, transform):
        """
        Helper function, called by _get_files() to process a plotter_transform.

        The 'transform' argument should be the json object corresponding to a plotter_transform,
        i.e. transform['name'] or transform['class_name'] should be 'plotter_transform...'.
        """
        self._add_plot_group(transform['plots'])


    def _add_bonsai_dedisperser(self, transform):
        """
        Helper function, called by _get_files() to process a bonsai_dedisperser transform.

        The 'transform' argument should be the json object corresponding to a bonsai_dedisperser,
        i.e. transform['name'] or transform['class_name'] should be 'bonsai_dedisperser...'.
        """

        # Here, we establish whether we're using the regular plotter or old bonsai plotter (in which case, we don't need to 
        # worry about separating part of the outputs into different rows) or the new bonsai plotter

        if ('n_plot_groups' not in transform):
            self._add_plot_group(transform['plots'])
            return
        
        # transform['plots'] will be a list of length (nloops * nzoom)
        assert len(transform['plots']) % transform['n_plot_groups'] == 0
        nloops = transform['n_plot_groups']
        group_size = len(transform['plots']) // nloops

        for i in xrange(nloops):
            j = (nloops-i-1) * group_size
            k = (nloops-i) * group_size
            self._add_plot_group(transform['plots'][j:k])

    
    def _add_plot_group(self, tp):
        """
        Helper function, called by _get_files() to process a new plot group.
        The 'tp' argument should be the 'plots' member of the json_data (a length of list nzoom).
        """
        
        # Start a new list for a new transform
        ftransform_group = []
        ttransform_group = []
        for zoom_level in tp:
            # This iterates over each zoom level (plot group) for a particular plotter transform (list of dictionaries)
            fzoom_group = []
            tzoom_group = []
            group_it0 = zoom_level['it0']
            for file_info in zoom_level['files'][0]:
                # We can finally access the file names :)
                name = file_info['filename']
                time = (group_it0 + file_info['it0']) * self.s_per_sample + self.t0   # start time of the plot in seconds
                fzoom_group.append(name)
                tzoom_group.append(time)
            ftransform_group.append(fzoom_group)
            ttransform_group.append(tzoom_group)
        # The plotter_transform defines zoom_level 0 to be most-zoomed-in, and zoom_level (N-1) to be
        # most-zoomed-out. The web viewer uses the opposite convention, so we reverse the order here.
        ftransform_group.reverse()
        ttransform_group.reverse()
        self.fnames.append(ftransform_group)
        self.ftimes.append(ttransform_group)


    def _check_zoom(self):
        """Back in the dark days without ch_frb_rfi, you could force the bonsai plotter and the transform plotter
        to produce a different number of zoom levels! This doesn't make sense, so the web_viewer won't display
        these runs."""
        if self.fnames is None:
            return False
        a = len(self.fnames[0])
        for element in self.fnames:
            if len(element) != a:
                return False
        return True


    def check_set(self, zoom, index):
        """Checks whether a link should be added at the top of the page to the next set of images in the series."""

        # For whatever reason, there are differing number of plots for
        # different transforms of the same zoom. This only returns false
        # if there are absolutely no images left (i.e. it will return true
        # if there is only one image available at a particular zoom because
        # one transform happened to output more than the rest). This means
        # we need to check again when we are displaying each individual
        # image whether it exists.

        if zoom >= self.max_zoom or zoom < self.min_zoom or index < self.min_index or index >= max([element[zoom] for element in self.max_index]):
            return False
        return True


    def check_image(self, transform, zoom, index):
        """Checks whether a particular image is available (because some transforms seem to produce more plots than others)"""

        if zoom >= self.max_zoom or zoom < self.min_zoom or index < self.min_index or index >= self.max_index[transform][zoom]:
            return False
        return True


    def __str__(self):
        s = ''
        if self.fnames is None:
            return 'None\n'
        tf_counter = 0
        for tf_group in self.fnames:
            s+= '* TRANSFORM ' + str(tf_counter) + ' *\n'
            tf_counter += 1
            zoom_counter = 0
            for zoom_group in tf_group:
                s += 'ZOOM ' + str(zoom_counter) + '\n'
                zoom_counter += 1
                for file in zoom_group:
                    s += str(file) + ' '
                s += '\n'
            s += '\n\n'
        return s


class Crawler():
    """
    Searches the two top directories pointed to by plots (assumed to be users -> pipeline runs). 
    Parser() is called for each pipeline run, creating a dictionary with information about 
    each pipeline run for each user. 

    Separate class here because I thought it might be nice for it to get other interesting metadata
    at some point. Could just be added to Parser if not. 

    Internal state:
        self.pipeline_dir[user][rundir] -> Parser
    """

    def __init__(self, path='static/plots'):
        self.path = path  # web_viewer root directory

        # Initially empty!  To reduce startup time, we no longer enumerate all
        # rundirs on startup.  Instead, new Parsers are added to the dictionary
        # as they get "discovered" with web requests.

        self.pipeline_dir = { }


    def get_users(self):
        """
        Called by web viewer toplevel page, to get a list of all users with pipeline runs.
        Since new entries may have been added since the last call to get_users(), we listdir()
        the web_viewer root directory, rather than relying on cached results.
        """
        return [ d for d in listdir(self.path) if isdir(join(self.path,d)) ]


    def get_unsorted_runs(self, user):
        """Similar to get_users(): returns a list of all run directories for a given user."""

        d = join(self.path, user)
        return [ r for r in listdir(d) if exists(join(d,r,'rf_pipeline_0.json')) ]


    def get_sorted_runs(self, user):
        """
        Sort runs by prefix, then by date/time.

        The return value is built out of lists and tuples as follows:
            [ (prefix1, [run1, run2, run3, ...]), 
              (prefix2, [...]),
             ... ]
        """

        runs = dict()

        for run in self.get_unsorted_runs(user):
            prefix = run[:-18]
            if prefix not in runs:
                # We need to add a new key
                runs[prefix] = [run]
            else:
                # Add to existing list
                runs[prefix].append(run)

        return [ (prefix, sorted(runs[prefix])) for prefix in sorted(runs.keys()) ]


    def run_exists(self, user, run):
        return exists(join(self.path, user, run, 'rf_pipeline_0.json'))


    def get_parser(self, user, run):
        """Returns Parser object, or None on error."""

        # Check for cached Parser.
        if self.pipeline_dir.has_key(user) and self.pipeline_dir[user].has_key(run):
            return self.pipeline_dir[user][run]

        # Cached Parser does not exist.  Does the run exist on disk?
        if not self.run_exists(user, run):
            return None

        try:
            p = Parser(join(self.path, user, run))
        except:
            return None

        # Add Parser to cache
        self.pipeline_dir.setdefault(user, {})
        self.pipeline_dir[user][run] = p
        
        return p


    def __str__(self):
        s = ""
        # For not writing out a ridiculous amount of information when trying to debug
        for user in self.pipeline_dir:
            if user == "mburhanpurkar":
                s += '*' * 160 + '\n'
                s += "USER: %s\n" % user
                for run in self.pipeline_dir[user]:
                    s += '-' * 160 + '\n'
                    s +=  "RUN: %s\n" % run
                    s += self.pipeline_dir[user][run].__str__()
        return s

# where we will search for users/runs/plots
path = environ.get('WEB_VIEWER_ROOT', 'static/plots')
print 'web_viewer root: %s' % path

master_directories = Crawler(path)     # dirs contains a dictionary in the form {'user1': {'run1': Parser1, 'run2': Parser2, ...}, ...}


####################################################################################################


@app.route("/")
def index():
    """Home page! Links to each of the users' pipeline runs."""

    display = '<h3>Users</h3>'

    for user in master_directories.get_users():
        display += '<li><a href="%s">%s</a>\n' % (url_for('runs', user=user), user)

    display += '<p><a href="https://github.com/mburhanpurkar/web_viewer">Instructions / Help / Documentation</a></p>'
    return display


@app.route("/<string:user>/runs")
def runs(user):
    """Displays links to the pipeline runs for a particular user."""

    display = '<h3>%s\'s pipeline runs</h3>\n' % user
    display += '<p>[&nbsp;&nbsp;&nbsp;<a href="%s">Back to List of Users</a>&nbsp;&nbsp;&nbsp;]\n' % url_for('index')
    display += '<p><ul>\n'

    for (prefix, runs) in master_directories.get_sorted_runs(user):
        run_urls = [ ]
        for run in runs:
            u = url_for('show_tiles', user=user, run=run, zoom=0, index1=0, index2=3)
            u = '<a href="%s">%s</a>' % (u, run[-17:])
            run_urls.append(u)

        display += '<p> <li> <b>%s</b><br>' % prefix

        n = len(run_urls)
        for i in xrange(0,n,6):
            j = min(i+6,n)
            display += (' | '.join(run_urls[i:j]))
            display += '<br>'

    return display


@app.route("/<string:user>/<string:run>/get_tile/<string:fname>")
def get_tile(user, run, fname):
    """Returns single image tile."""
    dirname = '%s/%s/%s' % (path, user, run)
    return send_from_directory(dirname, fname)


@app.route("/<string:user>/<string:run>/show_tiles/<int:zoom>/<int:index1>/<int:index2>")
def show_tiles(user, run, zoom, index1, index2):
    """Tiled image viewer! Shows all of the plots produced from a pipeline run at different zooms 
    across varying time intervals. The range of pictures shown can be changed to any values in 
    the url (index1 is the index of the first image shown and index2 is the index of the last 
    and defaults are set to 0 and 4 for the link accessed from the home page). The numbers displayed
    are the time in seconds at the start of the plot."""

    p = master_directories.get_parser(user, run)

    if p is None:
        return "The run was not found, or its json file could not be parsed."

    if p.fnames is None:
        return 'No files found.'

    if p.max_index is None:
        s = 'The number of zoom levels produced by the pipeline plotter was unequal to the number of plots produced ' \
            'by the bonsai plotter. This pipeline run cannot be displayed.'
        return s

    display = '<h3>Displaying Plots %d-%d at Zoom %d</h3>' % (index1, index2, (p.max_zoom - zoom - 1))  # account for resversal of zoom order in plotter
    display += '<table cellspacing="0" cellpadding="0">'

    for transform in reversed(range(len(p.fnames))):    # reversed to show triggers first
        display += '<tr>'
        # First, add plot times 
        for index in range(index1, index2 + 1):
            if p.check_image(transform, zoom, index):
                display += '<td>%s</td>' % p.ftimes[transform][zoom][index]
        display += '</tr>'
        # Now, add the images
        for index in range(index1, index2 + 1):
            if p.check_image(transform, zoom, index):
                display += '<td><img src="%s"></td>' % url_for('get_tile', user=user, run=run, fname=p.fnames[transform][zoom][index])
        display += '</tr><tr><td>&nbsp;</td></tr>'

    # Links to user and user/run pages
    display += '<p><center>[&nbsp;&nbsp;&nbsp;<a href="%s">Back to Users List</a>&nbsp;&nbsp;&nbsp;<a href="%s">Back to Your Runs</a>&nbsp;&nbsp;&nbsp;<a href="%s">' \
               'Show Triggers</a>&nbsp;&nbsp;&nbsp;<a href="%s">Show Last Transform</a>&nbsp;&nbsp;&nbsp;]</center></p>' \
               % (url_for('index'), url_for('runs', user=user), url_for('show_triggers', user=user, run=run, zoom=0), 
                  url_for('show_last_transform', user=user, run=run, zoom=zoom))

    # Plots to be linked
    display += '<p> <center> [&nbsp;&nbsp;&nbsp;'

    if p.check_set(zoom, index1 - 1):
        display += '<a href="%s">%s</a>&nbsp;&nbsp;&nbsp;' % ((url_for('show_tiles',
                    user=user, run=run, zoom=zoom, index1=index1 - 1, index2=index2 - 1)), 'Prev Time')
    else:
        display += 'Prev Time&nbsp;&nbsp;&nbsp;'

    if p.check_set(zoom, index1 + 1):
        display += '<a href="%s">%s</a>&nbsp;&nbsp;&nbsp;' % ((url_for('show_tiles',
                    user=user, run=run, zoom=zoom, index1=index1 + 1, index2=index2 + 1)), 'Next Time')
    else:
        display += 'Next Time&nbsp;&nbsp;&nbsp;'

    if p.check_set(zoom, index1 - (index2 - index1)):
        display += '<a href="%s">%s</a>&nbsp;&nbsp;&nbsp;' % ((url_for('show_tiles',
                    user=user, run=run, zoom=zoom, index1=index1 - (index2 - index1), index2=index2 - (index2 - index1))), 'Jump Back')
    else:
        display += 'Jump Back&nbsp;&nbsp;&nbsp;'

    if p.check_set(zoom, index1 + (index2 - index1)):
        display += '<a href="%s">%s</a>&nbsp;&nbsp;&nbsp;' % ((url_for('show_tiles',
                    user=user, run=run, zoom=zoom, index1=index1 + (index2 - index1), index2=index2 + (index2 - index1))), 'Jump Forward')
    else:
        display += 'Jump Forward&nbsp;&nbsp;&nbsp;'

    # For making the zooming preserve column number
    if (index2 - index1) % 2 == 0:
        new_index1 = index1 * 2 + (index2 - index1) / 2
        new_index2 = index2 * 2 - (index2 - index1) / 2
    else:
        new_index1 = index1 * 2 + ceil(index2 - index1) / 2 + 1
        new_index2 = index2 * 2 - ceil(index2 - index1) / 2 + 1

    if p.check_set(zoom + 1, index1 * 2):
        display += '<a href="%s">%s</a>&nbsp;&nbsp;&nbsp;' % ((url_for('show_tiles',
                    user=user, run=run, zoom=zoom + 1, index1=int(new_index1), index2=int(new_index2))), 'Zoom In')
    else:
        display += 'Zoom In&nbsp;&nbsp;&nbsp;'

    # More column preservation
    if (index2 - index1) % 2 == 0:
        new_index1 = (index1 - (index2 - index1) / 2) / 2
        new_index2 = (index2 + (index2 - index1) / 2) / 2
    else:
        new_index1 = (index1 - ceil((index2 - index1) / 2)) / 2
        new_index2 = (index2 + (ceil((index2 - index1) / 2) + 1)) / 2

    if p.check_set(zoom - 1, index1 // 2):
        display += '<a href="%s">%s</a>&nbsp;&nbsp;&nbsp;' % ((url_for('show_tiles',
                    user=user, run=run, zoom=zoom - 1, index1=int(new_index1), index2=int(new_index2))), 'Zoom Out')
    else:
        display += 'Zoom Out&nbsp;&nbsp;&nbsp;'
    display += ']</p> </center>'
    return display


@app.route("/<string:user>/<string:run>/show_last_transform/<int:zoom>")
def show_last_transform(user, run, zoom):
    """Displays the plots for the last transform at a given zoom horizontally. The zoom level can be changed by 
    changing the value in the url. Currently just indexes the second last value in fnames."""

    p = master_directories.get_parser(user, run)

    if p is None:
        return "The run was not found, or its json file could not be parsed."

    if p.fnames is None:
        return 'No files found.'

    if p.max_index is None:
        s = 'The number of zoom levels produced by the pipeline plotter was unequal to the number of plots produced ' \
            'by the bonsai plotter. This pipeline run cannot be displayed.'
        return s

    triggerList = p.fnames[-2]
    display = '<h3>Displaying Last Transform Plots at Zoom %s</h3>' % (p.max_zoom - zoom - 1)
    display += '<p><center>[&nbsp;&nbsp;&nbsp;<a href="%s">Back to Users List</a>&nbsp;&nbsp;&nbsp;<a href="%s">Back to Your Runs</a>' \
               '&nbsp;&nbsp;&nbsp;]</center></p>' % (url_for('index'), url_for('runs', user=user))
    display += '<table cellspacing="0" cellpadding="0"><tr>'

    last_row = 0
    current_row = 0

    for i, trigger in enumerate(triggerList[zoom]):
        temp = url_for('get_tile', user=user, run=run, fname=trigger)
        if i > 1 and i < p.max_index[-1][zoom] - 2:
            temp_link = url_for('show_tiles', user=user, run=run, zoom=zoom, index1=i - 2, index2=i + 1)
            display += '<td><a href="%s"><img src="%s"></a></td>' % (temp_link, temp)
        else:
            display += '<td><img src="%s"></td>' % temp
        current_row += 1
        if (current_row - last_row) == 5:
            last_row = current_row
            display += '</tr><tr><td>&nbsp;</td></tr><tr>'
    display += '</tr></table>'
    return display


@app.route("/<string:user>/<string:run>/show_triggers/<int:zoom>")
def show_triggers(user, run, zoom):
    """Displays all trigger plots at a given zoom horizontally. The zoom level can be changed by changing the value in the url. 
    Currently just indexes the last value in fnames."""

    p = master_directories.get_parser(user, run)

    if p is None:
        return "The run was not found, or its json file could not be parsed."

    if p.fnames is None:
        return 'No files found.'

    if p.max_index is None:
        s = 'The number of zoom levels produced by the pipeline plotter was unequal to the number of plots produced ' \
            'by the bonsai plotter. This pipeline run cannot be displayed.'
        return s

    zoom = int(zoom)

    triggerList = p.fnames[-1]
    display = '<h3>Displaying Trigger Plots at Zoom %s</h3>' % (p.max_zoom - zoom - 1)
    display += '<p><center>[&nbsp;&nbsp;&nbsp;<a href="%s">Back to Users List</a>&nbsp;&nbsp;&nbsp;<a href="%s">Back to Your Runs</a>' \
               '&nbsp;&nbsp;&nbsp;]</center></p>' % (url_for('index'), url_for('runs', user=user))
    display += '<table cellspacing="0" cellpadding="0"><tr>'

    last_row = 0
    current_row = 0

    for i, trigger in enumerate(triggerList[zoom]):
        temp = url_for('get_tile', user=user, run=run, fname=trigger)
        if i > 1 and i < p.max_index[-1][zoom] - 2:
            temp_link = url_for('show_tiles', user=user, run=run, zoom=zoom, index1=i - 2, index2=i + 1)
            display += '<td><a href="%s"><img src="%s"></a></td>' % (temp_link, temp)
        else:
            display += '<td><img src="%s"></td>' % temp
        current_row += 1
        if (current_row - last_row) == 5:
            last_row = current_row
            display += '</tr><tr><td>&nbsp;</td></tr><tr>'
    display += '</tr></table>'
    return display
