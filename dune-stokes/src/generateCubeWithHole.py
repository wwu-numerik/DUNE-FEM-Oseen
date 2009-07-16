#!/usr/bin/python

import sys
import math
import os
import time

## global defines

# about the outer rectangle
rectangle_length_x = 2.0
rectangle_length_y = 1.0

# about the inner ellipse
ellipse_radius_x = 0.45
ellipse_radius_y = 0.25
ellipse_center_x = 0.5 #rectangle_length_x / 2.0;
ellipse_center_y = rectangle_length_y / 2.0;

# about the number of points to approximate
number_of_points_per_quarter = 4

# about the boundary ids
id_of_ellipse_faces = 2
id_of_bottom_rectangle_faces = 3
id_of_right_rectangle_faces = 4
id_of_top_rectangle_faces = 5
id_of_left_rectangle_faces = 6

# about the files to be written
triangle_filename = 'cube_with_hole_in_2d.poly'

## done with global defines


## function definitions

# calculate points on the inner circle,
# beginning bottom left
# returns list of points in 2d
def generate_ellipse( center_x, center_y, radius_x, radius_y, id, n ):
	ellipse = []
	for i in range( 0, n ) : # means i = 0, ..., n-1:
		arc = math.pi * ( ( 5.0 / 4.0 ) + ( ( 2.0 * i ) / n ) )
		sinus_arc = math.sin( arc )
		cosinus_arc = math.cos( arc )
		x1 = center_x + ( radius_x * cosinus_arc )
		x2 = center_y + ( radius_y * sinus_arc )
		point = [ x1, x2 ]
		ellipse.append( [ point, id] )
	return ellipse
	
# calculate points on outer rectangle
# beginning bottom left
# rather trivial
# returns list of points in 2d
# [ [ point_x, point_y ], id_of_face_starting_clockwise_at_this_point ]
def generate_rectangle( length_x, length_y, ids ) :
	rectangle = []
	rectangle.append( [ [ 0.0, 0.0 ], ids[0] ] )
	rectangle.append( [ [ length_x, 0.0 ], ids[1] ] )
	rectangle.append( [ [ length_x, length_y ], ids[2] ] )
	rectangle.append( [ [ 0.0, length_y ], ids[3] ] )
	return rectangle

# write to trianlge
def write_to_triangle( ellipse, rectangle, hole, triangle_filename ) :
	file = open( triangle_filename, 'w' )
	# header
	file.write( '# #########################################################\n' )
	file.write( '# %s\n' %( triangle_filename ) )
	file.write( '# written by generateCubeWithHole.py on %s\n' %( time.strftime('%Y/%m/%d %H:%M:%S') ) )
	total_number_of_vertices = len( ellipse ) + len( rectangle )
	total_number_of_faces = len( ellipse ) + len( rectangle )
	file.write( '# we have \t%i vertices\n' %( total_number_of_vertices ) )
	file.write( '#\t\t%i faces\n' %( total_number_of_faces ) )
	file.write( '#\t\t%i hole\n' %(1) )
	file.write( '# #########################################################\n' )
	file.write( '#\n' )
	# write the vertices
	file.write( '# these are the %i vertices\n' %( total_number_of_vertices ) )
	file.write( '%i 2 0 0\n' %( total_number_of_vertices ) )
	vertex_number = 0
	# of the inner ellipse
	for point_with_id_on_ellipse in ellipse :
		point_on_ellipse = point_with_id_on_ellipse[ 0 ]
		x = point_on_ellipse[ 0 ]
		y = point_on_ellipse[ 1 ]
		file.write( '\t%i\t%f\t%f\n' %( vertex_number, x, y ) )
		vertex_number += 1
	# and the outer rectangle
	for point_with_id_on_rectangle in rectangle :
		point_on_rectangle = point_with_id_on_rectangle[ 0 ]
		x = point_on_rectangle[ 0 ]
		y = point_on_rectangle[ 1 ]
		file.write( '\t%i\t%f\t%f\n' %( vertex_number, x, y ) )
		vertex_number += 1
	file.write( '# these were the %i vertices\n' %( total_number_of_vertices ) )
	file.write( '#\n' )
	# write the faces with boundary id
	file.write( '# these are the %i faces\n' %( total_number_of_faces ) )
	file.write( '%i 1\n' %(total_number_of_faces) )
	face_number = 0
	# of the inner ellipse	
	for point_with_id_on_ellipse in ellipse :
		id = point_with_id_on_ellipse[ 1 ]
		file.write( '\t%i\t%i\t%i\t%i\n' %( face_number, face_number % len( ellipse ), ( face_number +1 ) % len( ellipse ), id ) )
		face_number += 1
	# of the outer rectangle
	for point_with_id_on_rectangle in rectangle :
		id = point_with_id_on_rectangle[ 1 ]
		file.write( '\t%i\t%i\t%i\t%i\n' %( face_number, face_number % len( rectangle ) + len( ellipse ) , ( face_number +1 ) % len( rectangle ) + len( ellipse ), id ) )
		face_number += 1
	file.write( '# these were the %i faces\n' %( total_number_of_faces ) )
	file.write( '#\n' )
	# write the holes
	file.write( '# this is the hole\n' )
	file.write( '1\n' )
	file.write('\t0\t%f\t%f\n' %( hole[ 0 ], hole[ 1 ] ) )
	file.write( '# this was the hole\n' )
	file.write( '#\n' )
	# ending
	file.write( '# #########################################################\n' )
	file.write( '# end of file %s\n' %( triangle_filename ) )
	file.write( '# written by generateCubeWithHole.py on %s\n' %( time.strftime('%Y/%m/%d %H:%M:%S') ) )
	file.write( '# #########################################################\n' )

## done with function definitions


## main
number_of_points = 4 * number_of_points_per_quarter
points_on_ellipse = generate_ellipse( ellipse_center_x, ellipse_center_y, ellipse_radius_x, ellipse_radius_y, id_of_ellipse_faces, number_of_points )
points_on_rectangle = generate_rectangle( rectangle_length_x, rectangle_length_y, [ id_of_bottom_rectangle_faces, id_of_right_rectangle_faces, id_of_top_rectangle_faces, id_of_left_rectangle_faces ] )
write_to_triangle( points_on_ellipse, points_on_rectangle, [ ellipse_center_x, ellipse_center_y ], triangle_filename )
