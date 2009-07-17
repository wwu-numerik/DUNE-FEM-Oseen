#!/usr/bin/python

import sys
import math
import os
import time

## global defines

# about the rectangles
rectangle_length_x = 4.0
rectangle_length_y = 4.0
number_of_rectangles_x = 2
number_of_rectangles_y = 2

# about the ellipses
ellipse_radius_x = ( rectangle_length_x / 2.0 ) * 0.5
ellipse_radius_y = ( rectangle_length_y / 2.0 ) * 0.5
ellipse_center_x = rectangle_length_x / 2.0;
ellipse_center_y = rectangle_length_y / 2.0;

# about the number of points to approximate the ellipses
number_of_points_per_quarter = 1

# about the boundary ids
id_of_ellipse_faces = 2
id_of_bottom_rectangle_faces = 3
id_of_right_rectangle_faces = 4
id_of_top_rectangle_faces = 5
id_of_left_rectangle_faces = 6

# about the files to be written
triangle_filename = 'homogeneous_perforated_domain_2d.poly'

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
def write_to_triangle( ellipsoids, outer_rectangle, holes, triangle_filename ) :
	file = open( triangle_filename, 'w' )
	# header
	file.write( '# ########################################################################\n' )
	file.write( '# %s\n' %( triangle_filename ) )
	file.write( '# written by generateHomogeneousPerforatedDomain.py on %s\n' %( time.strftime('%Y/%m/%d %H:%M:%S') ) )
	# get the points
	total_number_of_vertices = 0
	total_number_of_faces = 0
	for ellipsoid in ellipsoids :
		total_number_of_vertices += len( ellipsoid )
		total_number_of_faces += len( ellipsoid )
	total_number_of_vertices += len( outer_rectangle )
	total_number_of_faces += len( outer_rectangle )
	file.write( '# we have \t%i vertices\n' %( total_number_of_vertices ) )
	file.write( '#\t\t%i faces\n' %( total_number_of_faces ) )
	file.write( '#\t\t%i hole(s)\n' %( len( holes ) ) )
	file.write( '# ########################################################################\n' )
	file.write( '#\n' )
	# write the vertices
	file.write( '# these are the %i vertices\n' %( total_number_of_vertices ) )
	file.write( '%i 2 0 0\n' %( total_number_of_vertices ) )
	vertex_number = 0
	# of the inner ellipsoids
	for ellipsoid in ellipsoids :
		for point_with_id_on_ellipse in ellipsoid :
			point_on_ellipse = point_with_id_on_ellipse[ 0 ]
			x = point_on_ellipse[ 0 ]
			y = point_on_ellipse[ 1 ]
			file.write( '\t%i\t%f\t%f\n' %( vertex_number, x, y ) )
			vertex_number += 1
	# and the outer rectangle
	for point_with_id_on_rectangle in outer_rectangle :
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
	# of the inner ellipsoids
	number_of_ellipsoids_written = 0
	for ellipsoid in ellipsoids :
		for point_with_id_on_ellipse in ellipsoid :
			id = point_with_id_on_ellipse[ 1 ]
			file.write( '\t%i\t%i\t%i\t%i\n' %( face_number, ( face_number % len( ellipsoid ) ) + number_of_ellipsoids_written, ( ( face_number +1 ) % len( ellipsoid ) ) + number_of_ellipsoids_written, id ) )
			face_number += 1
		number_of_ellipsoids_written += len( ellipsoid )
	# of the outer rectangle
	for point_with_id_on_rectangle in outer_rectangle :
		id = point_with_id_on_rectangle[ 1 ]
		file.write( '\t%i\t%i\t%i\t%i\n' %( face_number, face_number % len( outer_rectangle ) + ( total_number_of_vertices - len( outer_rectangle ) ) , ( face_number +1 ) % len( outer_rectangle ) + ( total_number_of_vertices - len( outer_rectangle ) ), id ) )
		face_number += 1
	file.write( '# these were the %i faces\n' %( total_number_of_faces ) )
	file.write( '#\n' )
	# write the holes
	file.write( '# this are the %i holes\n' %( len( holes ) ) )
	file.write( '%i\n' %( len( holes ) ) )
	for hole in holes :
		file.write('\t0\t%f\t%f\n' %( hole[ 0 ], hole[ 1 ] ) )
	file.write( '# these were the %i holes\n' %( len( holes ) ) )
	file.write( '#\n' )
	# ending
	file.write( '# ########################################################################\n' )
	file.write( '# end of file %s\n' %( triangle_filename ) )
	file.write( '# written by generateHomogeneousPerforatedDomain.py on %s\n' %( time.strftime('%Y/%m/%d %H:%M:%S') ) )
	file.write( '# ########################################################################\n' )

## done with function definitions


## main

# points per ellipse
number_of_points = 4 * number_of_points_per_quarter

# generate the ellipsoids
ellipsoids = [ ]
holes = [ ]
for i in range( 0, number_of_rectangles_x ) :
  for j in range( 0, number_of_rectangles_y ) :
    center_x = ( i * rectangle_length_x ) + ellipse_center_x
    center_y = ( j * rectangle_length_y ) + ellipse_center_y
    ellipse = generate_ellipse( center_x, center_y, ellipse_center_x, ellipse_center_y, id_of_ellipse_faces, number_of_points )
    ellipsoids.append( ellipse )
    hole = [ center_x, center_y ]
    holes.append( hole )

# generate the outer rectangle
outer_rectangle = generate_rectangle( rectangle_length_x * number_of_rectangles_x, rectangle_length_y * number_of_rectangles_y, [ id_of_bottom_rectangle_faces, id_of_right_rectangle_faces, id_of_top_rectangle_faces, id_of_left_rectangle_faces ] )

# write to triangle .poly file
write_to_triangle( ellipsoids, outer_rectangle, holes, triangle_filename )

## done with main