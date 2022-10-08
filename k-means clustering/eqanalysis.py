"""
eqanalysis.py: analyze and plot earthquake data
Authors: Brandon Do
Credits:
Python Software Foundation: 'Turtle graphics for Tk'
https://docs.python.org/2/library/turtle.html

Youtube User Thenewboston: 'Python Programming Tutorial - 31 - Classes vs Instance Variables'
https://www.youtube.com/watch?v=qSDiHI1kP98

Youtube User Thenewboston: 'Python Programming Tutorial - 30 - init'
https://www.youtube.com/watch?v=G8kS24CtfoI

CIS 210 assignment 6, Fall 2015
"""

# a set of modules that we need to use in the code below
import math
import random
import argparse
from data import *
import turtle
import sys

# constants for the k-means clustering algorithm
# if you change these in your experimentation, you will need to look at 
# all parts of the code that refer to them, as there is some dependence
# on them (such as number of colors used in plotting clusters)
#
# IF YOU DO CHANGE THEM, YOU MUST PUT THEM BACK TO THE ORIGINAL VALUES
# BEFORE SUBMITTING YOUR WORK!!!!!
NO_OF_CLUSTERS = 6
NO_OF_ITERATIONS = 7

class Data_centen:
    """
        Uses list of specific data & provided functions for mean, median
        Args:
            eq_dict: dictionary of lists, each contained list represents an EQ event
            data_list: list of specific data to calculate central tendency
        Outputs:
            Statistics in a list
        """
    def __init__(self, data_list):
        self.data_list = data_list

    def statistics(self):
        mean = data_mean(self.data_list)
        median = data_median(self.data_list)

        return (mean, median)

class Data_disp:
    """
      Uses list of specific data & provided functions for variance & standard devation
        Args:
            eq_dict: dictionary of lists, each contained list represents an EQ event
            data_list: list of specific data to calculate dispersion statistics
        Outputs:
            Statistics in a list
        """
    def __init__(self, data_list):
        self.data_list = data_list

    def stdev(self):
        variance = data_mean_variance(self.data_list)
        stdev = math.sqrt(variance[1])
        return variance[1], stdev


class Data_iso:
    """
    Isolates columns of data, e.g. for magnitude it isolates the third column in the file,
    extracts them, and puts them into a list.
    Args:
        eq_dict: dictionary of lists, each contained list represents an EQ event
        column_number: location of the column. e.g. magnitude == column (2)
    Outputs:
        A list of the set of specific data, sorted.
    """
    def __init__(self, column_number, eq_dict):
        self.data_type = column_number
        self.eq_dict = eq_dict

    def isolator(self):
        iso_data = []
        for key in self.eq_dict:
            item = self.eq_dict[key]
            iso_data.append(item[self.data_type])
        iso_data.sort()
        return iso_data

    def xy_isolator(self):
        iso_data = []
        for key in self.eq_dict:
            item = self.eq_dict[key]
            xy_item = xy_calculate(item[0], item[1])
            iso_data.append([xy_item, item[self.data_type]])
        iso_data.sort()
        return iso_data


def euclid_distance(point1, point2):
    """
    computes the euclidean distance between two points
    Args:
        point1: list of floats, index 0 is longitude, index 1 is latitude
        point2: list of floats, index 0 is longitude, index 1 is latitude
    Returns:
        float, sqrt((x1-x2)**2 + (y1-y2)**2)
    """

    total = 0
    for index in range(2):
        diff = point1[index] - point2[index]
        total += diff * diff

    return math.sqrt(total)

def create_centroids(k, datadict):
    """
    randomly selects 'k' points from 'datadict' as the starting
        centroids for the k-means clustering algorithm
    Args:
        k: int, number of clusters desired
        datadict: list of lists, each contained list represents an EQ event
    Returns:
        list of lists, each contained list is an event to act as the centroid
    """
    centroids = []
    count = 0
    centroid_keys = []

    while count < k:
        rkey = random.randint(1, len(datadict))
        if rkey not in centroid_keys:
            centroids.append(datadict[rkey])
            centroid_keys.append(rkey)
            count += 1

    return centroids

def create_clusters(k, centroids, datadict, iterations):
    """
    k-means clustering algorithm - implementation taken from page 249 of
        ranum and miller text, with some modifications
    Args:
        k: integer, number of clusters
        centroids: list of events, each event is the centroid of its cluster
        datadict: dictionary of all EQ events
        iterations: int, number of clustering iterations to perform
    Returns:
        list of lists: each contained list is the set of indices into 'datadict'
           for events that belong to that cluster
    """
    for iteration in range(iterations):
        #print("****Iteration", iteration, "****")
        clusters = []
        for i in range(k):
            clusters.append([])

        for key in datadict:
            distances = []
            for cl_index in range(k):
                dist = euclid_distance(datadict[key], centroids[cl_index])
                distances.append(dist)
            min_dist = min(distances)
            index = distances.index(min_dist)
            clusters[index].append(key)

        dimensions = 2
        for cl_index in range(k):
            sums = [0]*dimensions
            for key in clusters[cl_index]:
                data_points = datadict[key]
                for ind in range(2):
                    sums[ind] = sums[ind] + data_points[ind]
            for ind in range(len(sums)):
                cl_len = len(clusters[cl_index])
                if cl_len != 0:
                    sums[ind] /= cl_len
            centroids[cl_index] = sums

        #for c in clusters:
            #print("CLUSTER")
            #for key in c:
                #print(datadict[key], end=" ")
            #print()

    return clusters

def read_file(filename):
    """
    read the EQ events from the csv file, 'filename'; any lines starting with
        # are skipped; the longitude, latitude, magnitude, and depth (in miles)
        is extracted from each event record, and stored as a list against its
        record number in a dictionary
    Args:
        filename: string, name of a CSV file containing the EQ data
    Returns:
        dictionary, indexed by integers, each value is a list of floats
            representing an EQ event
    """
    dict = {}
    key = 0

    fd = open(filename, "r")
    for line in fd:
        if line[0] == '#':
            continue		# causes the loop to grab another line
        key += 1
        values = line.rstrip('\n').split(',')
        lat = float(values[7])
        lon = float(values[8])
        mag = float(values[1])
        dep = float(values[10])
        dict[key] = [lon, lat, mag, dep]
    fd.close()
    return dict

# global data for map - if we had ;earmed about classes yet, this would have
# been hidden in a class instance, and the plot_*() functions would be methods
# on that class instance.  for now, these are global variables, and the
# plot functions access them

eq_turtle = None
eq_win = None
# these are the longitudes and latitudes for the Pacific NorthWest map that
# I have provided to you; do not change them!
left_lon = -128.608689
right_lon = -114.084764
top_lat = 51.248522
bot_lat = 38.584004
lon_diff = 0
lat_diff = 0
size_x = 0
size_y = 0
left_x = 0
bot_y = 0

def prepare_turtle():
    """
    Prepares the turtle and the window to plot magnitudes, depths, or clusters
    Args:
        None
    Outputs:
        creates turtle, sets window size, defines remainder of global
        data needed for plot_routines
    """
    global eq_turtle, eq_win
    global left_lon, right_lon, top_lat, bot_lat
    global lon_diff, lat_diff
    global size_x, size_y, left_x, bot_y

    eq_turtle = turtle.Turtle()
    eq_turtle.speed(10)
    eq_win = turtle.Screen()
    eq_win.screensize(655,808)	# number of pixels in the map I have provided
    lon_diff = right_lon - left_lon
    lat_diff = top_lat - bot_lat
    size_x = eq_win.screensize()[0]
    left_x = -size_x/2
    size_y = eq_win.screensize()[1]
    bot_y = -size_y/2
    eq_win.bgpic("PacificNW.gif")	# the map I have provided
    eq_turtle.hideturtle()
    eq_turtle.up()

def xy_calculate(lon, lat):
    """
    compute (x, y) given lon[gitude] and lat[itude]
    Args:
        lon: float, longitude value for point on map
        lat: float, latitude value for point on map
    Returns:
        tuple, corresponding pixel x and y values for use in turtle methods
    """
    global left_lon, right_lon, top_lat, bot_lat
    global lon_diff, lat_diff
    global size_x, size_y, left_x, bot_y

    x = left_x + (lon - left_lon) / lon_diff * size_x
    y = bot_y + (lat - bot_lat) / lat_diff * size_y
    return (x, y)

def plot_clusters(eq_clusters, eq_dict):
    """
    plot the clusters - use turtle.dot() at the appropriate location on the
        map for each event; use a different color for the events in each
        cluster - e.g. for cluster 0, use 'red', for 1, use 'violet' ...
    Args:
        eq_clusters: list of lists, each contained list has the indices for
                     events in that cluster in eq_dict
        eq_dict: dictionary of lists, each contained list represents an EQ event
    Outputs:
        plots all events in a particular cluster as dots on the map
    """

    COLORS = {
            0:'green', 1:'red', 2:'blue', 3:'cyan', 4:'violet',
            5:'purple', 6:'brown',7:'yellow', 8:'navy', 9:'light green',}

    global eq_turtle
    count = 0
    final_dict = {}
    for cluster in eq_clusters:
        cluster_list = []
        for number in cluster:
            item = eq_dict[number]
            xy_coord = xy_calculate(item[0], item[1])
            cluster_list.append(xy_coord)
        final_dict[count] = cluster_list
        count += 1

    for key in final_dict:
        item = final_dict[key]
        for xycoord in item:
            eq_turtle.goto(xycoord)
            eq_turtle.dot(7.5, COLORS[key])

def bin_value(value, bounds):
    """
    'bounds' defines a set of bins; this function returns the index of the
        first bin that contains 'value'
    Args:
        value: float, value to place in bin
        bounds: list of floats, bounds[i] is the top value of the bin
                code assumes that bounds is an increasing set of values
    Returns:
        integer, index of smallest value of bounds[] that is >= value
            if value > bounds[-1], returns len(bounds)
    """
    for i in range(len(bounds)):
        if value <= bounds[i]:
            return i
    return len(bounds)

def plot_magnitudes(eq_dict):
    """
    plot the magnitudes - use turtle.dot() at the appropriate location on the
        map for each event; use a different color and size for magnitude
        equivalence classes - e.g. if magnitude of event is <=1, use small dot
        that is 'violet', if between 1 and 2, use slightly larger dot that is
        'blue', ..., if between 9-10, use very large dot that is 'red'
    Args:
        eq_dict: dictionary of lists, each contained list represents an EQ event
    Outputs:
        plots magnitude of all events as dots on the map
    """

    global eq_turtle

    extraction_class = Data_iso(2, eq_dict)
    xypoint_list = extraction_class.xy_isolator()

    for point in xypoint_list:
        eq_turtle.goto(point[0])
        if point[1] <= 1.0:
            eq_turtle.dot(7.5, 'violet')
        if 1.0 < point[1] <= 2.0:
            eq_turtle.dot(15, 'blue')
        if point[1] > 9.0:
            eq_turtle.dot(22.5, 'red')

def plot_depths(eq_dict):
    """
    plot the depths - use turtle.dot() at the appropriate location on the
        map for each event; use a different color and size for depth
        equivalence classes - e.g. if depth of event is <=1 mile, use a large
        dot that is 'red', if between 1 and 5, use slightly smaller dot that is
        'orange', ..., if between 50-100, use a small dot that is 'violet'
    Args:
        eq_dict: dictionary of lists, each contained list represents an EQ event
    Outputs:
        plots depth of all events as dots on the map
    """
    global eq_turtle
    extraction_class = Data_iso(3, eq_dict)
    xypoint_list = extraction_class.xy_isolator()

    for point in xypoint_list:
        eq_turtle.goto(point[0])
        if 50 <= point[1] <= 100:
            eq_turtle.dot(5, 'cyan')
        if 1 <= point[1] <= 5:
            eq_turtle.dot(10, 'blue')
        if point[1] <= 1:
            eq_turtle.dot(15, 'green')


def analyze_depths(eq_dict):
    """
    Perform statistical analysis on the depth information in the dictionary
    Args:
        eq_dict: dictionary of lists, each contained list represents an EQ event
    Outputs:
        mean, median, and standard deviation of depth data
        frequency table for the depth data
    """

    depth_list = Data_iso(3, eq_dict)
    depth_centen = Data_centen(depth_list.isolator())
    depth_disp = Data_disp(depth_list.isolator())

    centen = depth_centen.statistics() #0 - mean, 1 - median
    disp = depth_disp.stdev()          #0 - variance, 1 - standard deviation

    frequency = frequency_list(depth_list.isolator())

    data_format(centen[0], centen[1], disp[1], 'depth', frequency, 'miles')

def data_format(mean, median, stdev, type_name, frequency, units):
    """
    Takes the central tendency statistics and prints in in a format
    that is viewable for the user. Used for 'analyze magnitude/depths'
    command.
    Args:
        mean, median, stdev: Results of inputting the data list into the given
                                data.py functions
        frequency: A frequency table from the 'frequency_list' function
        type_name: The name of the type of data, ie 'magnitude'
        units: Type of number to add to print, ie 'Miles', if none then units=''
    Outputs:
        No output, prints all the data for the user.
    """
    print('Analysis of {} data'.format (type_name))
    print('     Mean {} = {:.02} {}'.format(type_name, mean, units))
    print('     Median {} = {:.02} {}'.format(type_name, median, units))
    print('     Standard Deviation = {:.02f} {}'.format(stdev, units))
    print('ITEM    FREQUENCY')
    for item in frequency:
        print('  {}      {}'.format(item[0], item[1]))

def frequency_list(data_list):
    """
    Takes a list of data of one catagory (ie magnitude) and counts frequency
    of each occurance. Each unique point is stored in a list with its count.
    Used for 'analyze depths/magnitudes' command.
    Args:
        data_list: Data list of one category
    Outputs:
        List of lists which have two elements each, the unique key and frequency. ie [[zealot, 2], [tracer, 4]]
    """
    
    freq_list = []
    checked_list = []

    for item in data_list:
        if item not in checked_list:
            freq_list.append([item, data_list.count(item)])
            checked_list.append(item)

    return freq_list

def analyze_magnitudes(eq_dict):
    """
    Perform statistical analysis on the magnitude information in the dictionary
    Args:
        eq_dict: dictionary of lists, each contained list represents an EQ event
    Outputs:
        mean, median, and standard deviation of magnitude data
        frequency table for the magnitude data
    """

    magnitude_list = Data_iso(2, eq_dict)
    magnitude_centen = Data_centen(magnitude_list.isolator())
    magnitude_disp = Data_disp(magnitude_list.isolator())

    centen = magnitude_centen.statistics() #0 - mean, 1 - median
    disp = magnitude_disp.stdev()          #0 - variance, 1 - standard deviation

    frequency = frequency_list(magnitude_list.isolator())
    units = ''

    data_format(centen[0], centen[1], disp[1], 'Magnitude', frequency, units)

def analyze_clusters(eq_clusters, eq_dict):
    """
    Perform statistical analysis on the depth and magnitude information
        for each cluster
    Args:
        eq_clusters: list of lists, each contained list has the indices into
                     eq_dict for events in that cluster
        eq_dict: dictionary of lists, each contained list represents an EQ event
    Outputs:
        put into a dictionary, for each cluster:
            mean, median, and standard deviation of magnitude data
            mean, median, and standard deviation of depth data
    """

    counter = 0
    for cluster in eq_clusters:
        cluster_dict = {}

        for number in cluster:
            cluster_dict[number] = eq_dict[number]

        #Cluster Statistics for Magnitude
        magnitude_list = Data_iso(2, cluster_dict)
        magnitude_centen = Data_centen(magnitude_list.isolator())
        magnitude_disp = Data_disp(magnitude_list.isolator())
        m_centen = magnitude_centen.statistics()  # 0 - mean, 1 - median
        m_disp = magnitude_disp.stdev()  # 0 - variance, 1 - standard deviation

        #Cluster Statistics for Depth
        depth_list = Data_iso(3, cluster_dict)
        depth_centen = Data_centen(depth_list.isolator())
        depth_disp = Data_disp(depth_list.isolator())
        d_centen = depth_centen.statistics()  # 0 - mean, 1 - median
        d_disp = depth_disp.stdev()  # 0 - variance, 1 - standard deviation

        #Prints the data for the user, repeats for cluster.
        print('Analysis of cluster {}'.format(counter))
        print('     Analysis of magnitude data')
        print('         Mean magnitude = {:.1f}'.format(m_centen[0]))
        print('         Median magnitude = {:.1f}'.format(m_centen[1]))
        print('         Standard deviation = {:.02f}'.format(m_disp[1]))
        print('     Analysis of depth data')
        print('         Mean depth = {:.1f} miles'.format(d_centen[0]))
        print('         Median depth = {:.1f} miles'.format(d_centen[1]))
        print('         Standard deviation = {:.02f} miles'.format(d_disp[1]))

        counter += 1

def main():
    """
    Interaction if run from the command line.
    Usage:  python3 eqanalysis.py eq_data_file.csv command
    """
    parser = argparse.ArgumentParser(description="Earthquake event file stats")
    parser.add_argument('eq_file', type=str,
                 help='A csv file containing earthquake events, one per line.')
    parser.add_argument('command', type=str,
                 help='One of the following strings: plot analyze')
    parser.add_argument('what', type=str,
                 help='One of the following strings: clusters depths magnitudes')
    args = parser.parse_args()
    eq_file = args.eq_file
    cmd = args.command
    what = args.what
    if cmd != 'plot' and cmd != 'analyze':
        print('Illegal command: {}; must be "plot" or "analyze"'.format(cmd))
        sys.exit(1)
    if what != 'clusters' and what != 'magnitudes' and what != 'depths':
        print('Can only process clusters, magnitudes, or depths')
        sys.exit(1)
    eq_dict = read_file(eq_file)
    prepare_turtle()
    if what == 'clusters':
        eq_centroids = create_centroids(NO_OF_CLUSTERS, eq_dict)
        eq_clusters = create_clusters(NO_OF_CLUSTERS, eq_centroids, eq_dict, NO_OF_ITERATIONS)
    if cmd == 'plot':
        if what == 'clusters':
            plot_clusters(eq_clusters, eq_dict)
        elif what == 'magnitudes':
            plot_magnitudes(eq_dict)
        elif what == 'depths':
            plot_depths(eq_dict)
        print("ALL EVENTS HAVE BEEN PLOTTED")
        eq_win.exitonclick()
    else:
        if what == 'clusters':
            analyze_clusters(eq_clusters, eq_dict)
        elif what == 'magnitudes':
            analyze_magnitudes(eq_dict)
        elif what == 'depths':
            analyze_depths(eq_dict)

if __name__ == "__main__":
    main()