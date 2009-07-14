#!/usr/bin/python

import sys
import math
import os

## global defines

# about the outer rectangle
rectangle_length_x = 2.0
rectangle_length_y = 1.0

# about the inner ellipse
ellipse_radius_x = 0.2
ellipse_radius_y = 0.1
ellipse_center_x = rectangle_length_x / 2.0;
ellipse_center_y = rectangle_length_y / 2.0;

# about the number of points to approximate
number_of_points_per_quarter = 2

# about the boundary ids
id_of_circle_faces = 2
id_of_bottom_square_faces = 3
id_of_right_square_faces = 4
id_of_top_square_faces = 5
id_of_left_square_faces = 6

# about the fiels to be written
dgf_filename = 'cube_with_hole_in_2d.dgf'

# about triangle

## done with global defines


## function definitions

# calculate points on the inner circle,
# beginning bottom left
def generate_ellipse( center_x, center_y, radius_x, radius_y, n ):
	ellipse = []
	for i in range( 0, n ) : # means i = 0, ..., n-1:
		arc = math.pi * ( ( 5.0 / 4.0 ) + ( ( 2.0 * i ) / n ) )
		sinus_arc = math.sin( arc )
		cosinus_arc = math.cos( arc )
		x1 = center_x + ( radius_x * cosinus_arc )
		x2 = center_y + ( radius_y * sinus_arc )
		point = [ x1, x2 ]
		ellipse.append( point )
	return ellipse
	
# calculate points on outer rectangle
# beginning bottom left
# rather trivial
def generate_rectangle( length_x, length_y ) :
	rectangle = []
	rectangle.append( [ 0.0, 0.0 ] )
	rectangle.append( [ length_x, 0.0 ] )
	rectangle.append( [ 0.0, length_y ] )
	rectangle.append( [ length_x, length_y ] )
	return rectangle

# write 

## done with function definitions
	



# main
#read_command_line_arguments()
number_of_points = 4 * number_of_points_per_quarter
points_on_ellips = generate_ellipse( ellipse_center_x, ellipse_center_y, ellipse_radius_x, ellipse_radius_y, number_of_points )
points_on_rectangle = generate_rectangle( rectangle_length_x, rectangle_length_y )
print points_on_rectangle
#generate_inner_circle_point_list( number_of_points )
#generate_square_around_cirle_point_list( number_of_points )
#write_to_dgf_file( number_of_points )

## get command line arguments
#def read_command_line_arguments():
    #global circle_radius
    #global square_length
    #global number_of_points_per_side
    #if len( sys.argv ) == 4:
        #circle_radius = float( sys.argv[ 1 ] )
        #square_length = float( sys.argv[ 2 ] )
        #number_of_points_per_side = int( sys.argv[ 3 ] )
        #print 'circle_radius:', circle_radius
        #print 'square_length:', square_length
        #print 'number_of_points_per_side:', number_of_points_per_side
    #elif len( sys.argv ) == 2:
        #number_of_points_per_side = int( sys.argv[ 1 ] )
        #print 'circle_radius:', circle_radius
        #print 'square_length:', square_length
        #print 'number_of_points_per_side:', number_of_points_per_side
    #elif len( sys.argv ) == 1:
        #print 'circle_radius:', circle_radius
        #print 'circle_radius_x:', circle_radius_x
        #print 'circle_radius_y:', circle_radius_y
        #print 'square_length:', square_length
        #print 'number_of_points_per_side:', number_of_points_per_side
    #else:
        #print 'usage:'
        #print '\t', sys.argv[ 0 ], 'number_of_points_per_side'
        #print 'or\t', sys.argv[ 0 ], 'circle_radius square_length number_of_points_per_side'
		
## calculate points on the outer square,
#def generate_square_around_cirle_point_list( n ):
    #global points_on_square
	## bottom face first
    #for i in range( 0, n / 4 ):
        #shift = square_length_x / ( n / 4.0 )
        #point_on_square = [ i * shift, 0.0 ]
        #points_on_square.append( point_on_square )
    ## then right face
    #for i in range( n / 4, n / 2 ):
        #shift = square_length_y / ( n / 4.0 )
        #point_on_square = [ square_length_x, shift * ( i - ( n / 4.0 ) ) ]
        #points_on_square.append( point_on_square )
    ## then top face
    #for i in range( n / 2, 3 * ( n / 4 ) ):
        #shift = square_length_x / ( n / 4.0 )
        #point_on_square = [ square_length_x - shift * ( i - n / 2.0 ), square_length_y ]
        #points_on_square.append( point_on_square )
    ## then left face
    #for i in range( 3 * ( n / 4 ), n ):
        #shift = square_length_y / ( n / 4.0 )
        #point_on_square = [ 0.0, square_length_y - shift * ( i - ( 3.0 * n ) / 4.0 ) ]
        #points_on_square.append( point_on_square )

## write to dgf file
#def write_to_dgf_file( n ):
    #dgf_file = open( dgf_filename, 'w' )
    #dgf_file.write( 'DGF\n' )
    #dgf_file.write( '#\n' )
    #dgf_file.write( 'VERTEX\t\t\t# the vertices\n' )
    #number_of_vertices = 0
    #for point in points_on_square:
        #dgf_file.write( '%f %f\t# vertex %i,\t square\n' %( point[0], point[1], number_of_vertices ) )
        #number_of_vertices += 1
    #for point in points_on_circle:
        #dgf_file.write( '%f %f\t# vertex %i,\t circle\n' %( point[0], point[1], number_of_vertices ) )
        #number_of_vertices += 1
    #dgf_file.write( '#\n' )
    #dgf_file.write( 'CUBE\t\t\t# the cubes\n' )
    #number_of_cubes = 0
    #for i in range( 0, n ):
        #dgf_file.write( '%i %i %i %i\t\t\t# cube %i,\t from vertices %i, %i, %i and %i\n' %( i, ( i + 1 ) % n, n + i, ( n + i + 1 ) % n + n, number_of_cubes, i,  ( i + 1 ) % n, n + i, ( n + i + 1 ) % n + n ) )
        #number_of_cubes += 1
    #dgf_file.write( '#\n' )
    #dgf_file.write( 'BOUNDARYSEGMENTS\t\t# the boundary ids\n' )
    #vertex_number = 0
    #for i in range( 0, n / 4 ):
        #dgf_file.write( '%i %i %i\t\t\t# bottom face of the square\n' %( id_of_bottom_square_faces, vertex_number, ( vertex_number + 1 ) % n ) )
        #vertex_number += 1
    #for i in range( n / 4, n / 2 ):
        #dgf_file.write( '%i %i %i\t\t\t# right face of the square\n' %( id_of_right_square_faces, vertex_number, ( vertex_number + 1 ) % n ) )
        #vertex_number += 1
    #for i in range( n / 2, 3 * ( n / 4 ) ):
        #dgf_file.write( '%i %i %i\t\t\t# top face of the square\n' %( id_of_top_square_faces, vertex_number, ( vertex_number + 1 ) % n ) )
        #vertex_number += 1
    #for i in range( 3 * ( n / 4 ), n ):
        #dgf_file.write( '%i %i %i\t\t\t# left face of the square\n' %( id_of_left_square_faces, vertex_number, ( vertex_number + 1 ) % n ) )
        #vertex_number += 1
    #for i in range( n, 2 * n ):
        #dgf_file.write( '%i %i %i\t\t\t# of the face of the circle\n' %( id_of_circle_faces, i, ( i + 1 ) % n + n) )
    #dgf_file.write( '#\n' )
    #dgf_file.write( 'BOUNDARYDOMAIN\n' )
    #dgf_file.write( 'default 1\t\t# default value for all other faces (which should not exist)\n' )
    #dgf_file.write( '#\n' )
