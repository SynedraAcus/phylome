#! /usr/bin/env python3

from argparse import ArgumentParser
from collections import namedtuple
import re
import svgwrite


class Point:
    def __init__(self, x, y, a, b, c):
        self.x = x
        self.y = y
        self.a = a
        self.b = b
        self.c = c


def barycentric_to_euclidean(point, triangle):
    """
    Convert barycentric coordinates in the triangle to euclidean coords with
    the origin in the triangle's center
    :param bary:
    :return:
    """
    point.x = point.a * triangle[0].x + point.b * triangle[1].x + \
              point.c * triangle[2].x
    point.y = point.a * triangle[0].y + point.b * triangle[1].y + \
              point.c * triangle[2].y


def barycentric_polygon(points, color, triangle):
    """
    Return svgwrite-compatible Polygon object with given points
    :param points:
    :param color:
    :return:
    """
    for point in points:
        barycentric_to_euclidean(point, triangle)
    return svgwrite.shapes.Polygon(points=[(point.x, point.y) for point in points],
                                   fill=color, stroke='black')


def percentage_color(percentage):
    """
    Return color string in 'rrggbb' notation, where 0% is ffffff
    and 100% is 000000
    :param percentage:
    :return:
    """
    assert 0 <= percentage <= 100
    s = str(str(hex(round(255*(1-percentage/100)))))
    return '#'+s.split('x')[1]*3


parser = ArgumentParser('Process the IQ-tree likelihood mapping data')
parser.add_argument('-f', nargs='*', help='iqtree files')
args = parser.parse_args()

clusters = ('diatom', 'red', 'green', 'rest')
cutoff_line = '-'*144+'\n'
regex = re.compile('\(([\d.]+) *\)')
segment_values = []
count = 0
for iqfile in args.f:
    # Read the summary segment support from seven-segment table
    # It'd be right after the table, which begins and ends with cutoff line
    wait = 0
    for line in open(iqfile):
        if wait == 2:
            segment_values.append([float(x) for x in re.findall(regex, line)])
            count += 1
            break
        if line == cutoff_line:
            wait += 1

summary = [0 for x in range(7)]
for cluster_run in segment_values:
    for index, value in enumerate(cluster_run):
        summary[index] += value
summary = [x/count for x in summary]
print(summary)

# Three corners with euclidean and barycentric coordinates
triangle = [Point(500, 100, 1, 0, 0),
            Point(846, 700, 0, 1, 0),
            Point(154, 700, 0, 0, 1)]
drawing = svgwrite.Drawing(filename='Triangle.svg',
                           size=(1000, 1000))
# Background
drawing.add(svgwrite.shapes.Rect(insert=(0, 0), size=(1000, 1000),
                                 fill='white'))
# Code below could be written better, but whatever
# Segment 1
drawing.add(barycentric_polygon([Point(None, None, 1, 0, 0),
                                 Point(None, None, 3/4, 1/4, 0),
                                 Point(None, None, 2/3, 1/6, 1/6),
                                 Point(None, None, 3/4, 0, 1/4)],
                                 percentage_color(summary[0]), triangle))
text_point = Point(None, None, 5/6, 1/12, 1/12)
barycentric_to_euclidean(text_point, triangle)
drawing.add(svgwrite.text.Text(str(summary[0]) + ' %',
                               insert=(text_point.x, text_point.y),
                               font_size=20, text_anchor='middle'))
# Segment 2
drawing.add(barycentric_polygon([Point(None, None, 1/6, 2/3, 1/6),
                                 Point(None, None, 1/4, 3/4, 0),
                                 Point(None, None, 0, 1, 0),
                                 Point(None, None, 0, 3/4, 1/4)],
                                 percentage_color(summary[1]), triangle))
text_point = Point(None, None, 1/12, 5/6, 1/12)
barycentric_to_euclidean(text_point, triangle)
drawing.add(svgwrite.text.Text(str(summary[1]) + ' %',
                               insert=(text_point.x, text_point.y),
                               font_size=20, text_anchor='middle'))

# Segment 3
drawing.add(barycentric_polygon([Point(None, None, 1/4, 0, 3/4),
                                 Point(None, None, 1/6, 1/6, 2/3),
                                 Point(None, None, 0, 1/4, 3/4),
                                 Point(None, None, 0, 0, 1)],
                                percentage_color(summary[2]), triangle))
text_point = Point(None, None, 1/12, 1/12, 5/6)
barycentric_to_euclidean(text_point, triangle)
drawing.add(svgwrite.text.Text(str(summary[2]) + ' %',
                               insert=(text_point.x, text_point.y),
                               font_size=20, text_anchor='middle'))

# Segment 4
drawing.add(barycentric_polygon([Point(None, None, 3/4, 1/4, 0),
                                 Point(None, None, 1/4, 3/4, 0),
                                 Point(None, None, 1/6, 2/3, 1/6),
                                 Point(None, None, 2/3, 1/6, 1/6)],
                                percentage_color(summary[3]), triangle))
text_point = Point(None, None, 11/24, 11/24, 1/12)
barycentric_to_euclidean(text_point, triangle)
drawing.add(svgwrite.text.Text(str(summary[3]) + ' %',
                               insert=(text_point.x, text_point.y),
                               font_size=20, text_anchor='middle'))

# Segment 5
drawing.add(barycentric_polygon([Point(None, None, 1/6, 2/3, 1/6),
                                 Point(None, None, 0, 3/4, 1/4),
                                 Point(None, None, 0, 1/4, 3/4),
                                 Point(None, None, 1/6, 1/6, 2/3)],
                                percentage_color(summary[4]), triangle))
text_point = Point(None, None, 1/12, 11/24, 11/24)
barycentric_to_euclidean(text_point, triangle)
drawing.add(svgwrite.text.Text(str(summary[4]) + ' %',
                               insert=(text_point.x, text_point.y),
                               font_size=20, text_anchor='middle'))

# Segment 6
drawing.add(barycentric_polygon([Point(None, None, 3/4, 0, 1/4),
                                 Point(None, None, 2/3, 1/6, 1/6),
                                 Point(None, None, 1/6, 1/6, 2/3),
                                 Point(None, None, 1/4, 0, 3/4)],
                                percentage_color(summary[5]), triangle))
text_point = Point(None, None, 11/24, 1/12, 11/24)
barycentric_to_euclidean(text_point, triangle)
drawing.add(svgwrite.text.Text(str(summary[5]) + ' %',
                               insert=(text_point.x, text_point.y),
                               font_size=20, text_anchor='middle'))

# Segment 7
drawing.add(barycentric_polygon([Point(None, None, 2/3, 1/6, 1/6),
                                 Point(None, None, 1/6, 2/3, 1/6),
                                 Point(None, None, 1/6, 1/6, 2/3)],
                                percentage_color(summary[6]), triangle))
text_point = Point(None, None, 1/3, 1/3, 1/3)
barycentric_to_euclidean(text_point, triangle)
drawing.add(svgwrite.text.Text(str(summary[6]) + ' %',
                               insert=(text_point.x, text_point.y),
                               font_size=20, text_anchor='middle'))

drawing.save()
