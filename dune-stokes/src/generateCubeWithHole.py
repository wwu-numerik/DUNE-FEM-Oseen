#!/usr/bin/python

import sys
import math
import os

# global defines
circle_radius = 0.25
square_length = 1.0
number_of_points_per_side = 1
dgf_filename = 'cube_with_hole_in_2d.dgf'
id_of_circle_faces = 2
id_of_bottom_square_faces = 2
id_of_right_square_faces = 3
id_of_top_square_faces = 2
id_of_left_square_faces = 4

# get command line arguments
if len( sys.argv ) == 4:
	circle_radius = float( sys.argv[ 1 ] )
	square_length = float( sys.argv[ 2 ] )
	number_of_points_per_side = int( sys.argv[ 3 ] )
        print 'circle_radius:', circle_radius
        print 'square_length:', square_length
        print 'number_of_points_per_side:', number_of_points_per_side
elif len( sys.argv ) == 2:
	number_of_points_per_side = int( sys.argv[ 1 ] )
	print 'circle_radius:', circle_radius
	print 'square_length:', square_length
	print 'number_of_points_per_side:', number_of_points_per_side
elif len( sys.argv ) == 1:
        print 'circle_radius:', circle_radius
        print 'square_length:', square_length
        print 'number_of_points_per_side:', number_of_points_per_side
else:
	print 'usage:'
	print '\t', sys.argv[ 0 ], 'number_of_points_per_side'
	print 'or\t', sys.argv[ 0 ], 'circle_radius square_length number_of_points_per_side'

# some needed vaiables
n = 4 * number_of_points_per_side
points_on_circle = [  ]
points_on_square = [  ]

# calculate points on the inner circle,
# beginning in the top left bottom
for i in range( 0, n ): # means i = 0, ..., n-1
	arc = math.pi * ( ( 5.0 / 4.0 ) + ( ( 2.0 * i ) / n ) )
	sinus_arc = math.sin( arc )
	cosinus_arc = math.cos( arc )
	x1 = ( square_length / 2.0 ) + ( circle_radius * cosinus_arc )
	x2 = ( square_length / 2.0 ) + ( circle_radius * sinus_arc )
	point_on_circle = [ x1, x2 ]
	points_on_circle.append( point_on_circle )
	
# calculate points on the outer square,
# bottom face first
for i in range( 0, n / 4 ): 
	shift = square_length / ( n / 4.0 )
	point_on_square = [ i * shift, 0.0 ]
	points_on_square.append( point_on_square )
# then right face
for i in range( n / 4, n / 2 ): 
	shift = square_length / ( n / 4.0 )
	point_on_square = [ square_length, shift * ( i - ( n / 4.0 ) ) ]
	points_on_square.append( point_on_square )
# then top face
for i in range( n / 2, 3 * ( n / 4 ) ): 
	shift = square_length / ( n / 4.0 )
	point_on_square = [ square_length - shift * ( i - n / 2.0 ), square_length ]
	points_on_square.append( point_on_square )
# then left face
for i in range( 3 * ( n / 4 ), n ): 
	shift = square_length / ( n / 4.0 )
	point_on_square = [ 0.0, square_length - shift * ( i - ( 3.0 * n ) / 4.0 ) ]
	points_on_square.append( point_on_square )
	
# write to dgf file
dgf_file = open( dgf_filename, 'w' )
dgf_file.write( 'DGF\n' )
dgf_file.write( '#\n' )
dgf_file.write( 'VERTEX\t\t\t# the vertices\n' )
number_of_vertices = 0
for point in points_on_square:
	dgf_file.write( '%f %f\t# vertex %i,\t square\n' %( point[0], point[1], number_of_vertices ) )
	number_of_vertices += 1
for point in points_on_circle:
	dgf_file.write( '%f %f\t# vertex %i,\t circle\n' %( point[0], point[1], number_of_vertices ) )
	number_of_vertices += 1
dgf_file.write( '#\n' )
dgf_file.write( 'CUBE\t\t\t# the cubes\n' )
number_of_cubes = 0
for i in range( 0, n ):
	dgf_file.write( '%i %i %i %i\t\t\t# cube %i,\t from vertices %i, %i, %i and %i\n' %( i, ( i + 1 ) % n, n + i, ( n + i + 1 ) % n + n, number_of_cubes, i,  ( i + 1 ) % n, n + i, ( n + i + 1 ) % n + n ) )
	number_of_cubes += 1
dgf_file.write( '#\n' )
dgf_file.write( 'BOUNDARYSEGMENTS\t# the boundary ids\n' )
vertex_number = 0
for i in range( 0, n / 4 ): 
	dgf_file.write( '%i %i %i\t\t\t# bottom face of the square\n' %( id_of_bottom_square_faces, vertex_number, ( vertex_number + 1 ) % n ) )
	vertex_number += 1
for i in range( n / 4, n / 2 ): 
	dgf_file.write( '%i %i %i\t\t\t# right face of the square\n' %( id_of_right_square_faces, vertex_number, ( vertex_number + 1 ) % n ) )
	vertex_number += 1
for i in range( n / 2, 3 * ( n / 4 ) ): 
	dgf_file.write( '%i %i %i\t\t\t# top face of the square\n' %( id_of_top_square_faces, vertex_number, ( vertex_number + 1 ) % n ) )
	vertex_number += 1
for i in range( 3 * ( n / 4 ), n ): 
	dgf_file.write( '%i %i %i\t\t\t# left face of the square\n' %( id_of_left_square_faces, vertex_number, ( vertex_number + 1 ) % n ) )
	vertex_number += 1
for i in range( n, 2 * n ):
	dgf_file.write( '%i %i %i\t\t\t# of the face of the circle\n' %( id_of_circle_faces, i, ( i + 1 ) % n + n) )
dgf_file.write( '#\n' )
dgf_file.write( 'BOUNDARYDOMAIN\n' )
dgf_file.write( 'default 1\t\t# default value for all other faces (which should not exist)\n' )
dgf_file.write( '#\n' )
