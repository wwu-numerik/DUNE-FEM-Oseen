#!/usr/bin/python
# -*- coding: utf-8 -*-

from gridhelper import *
import math, copy

area 			= math.pi
tube_length		= 15
num_verts 		= 5
alpha			= math.radians(360. / (num_verts ) )
alpha_half		= alpha / 2.
num_midrings	= 3
area_one_tri = math.pi / float(num_verts)
L_x = math.sqrt( area_one_tri / ( math.sin( alpha_half ) * math.cos( alpha_half ) ) )

origin_L = Vector3(0,0,0)
points_L = PLCPointList( 3 )

bound_L = FanningSimplexList( points_L.appendVert( origin_L ), 2 )
L = Vector3( L_x, 0, 0 )
bound_L.addVertex( points_L.appendVert( L ) )
rot_mat = Matrix4.new_rotatez( alpha )
for i in range( 1, num_verts  ):
	L = rot_mat * L
	bound_L.addVertex(points_L.appendVert( L ))
bound_L.close()

grid = FullGrid( bound_L, 4 )

incr = tube_length / ( num_midrings + 1 )
for i in range( 1, num_midrings + 2 ):	
	origin_M = Vector3(0,0,incr*i )
	points_M = PLCPointList( 3 )
	bound_M = InbetweenRing( )
	M = Vector3( L_x, 0, incr*i )
	bound_M.addVertex( points_M.appendVert( M ) )
	for i in range( 1, num_verts  ):
		M = rot_mat  * M
		bound_M.addVertex(points_M.appendVert( M ))
	grid.connect(bound_M)
	
kipp_mat = Matrix4.new_rotatey( math.radians(15) )
origin_R = Vector3(0,0,tube_length )
points_R = PLCPointList( 3 )
bound_R = FanningSimplexList( points_R.appendVert( origin_R ), 3 )
R = Vector3( L_x, 0, tube_length )
#R = kipp_mat * R #skalier groesser R
bound_R.addVertex( points_R.appendVert( R ) )
for i in range( 1, num_verts  ):
	R = rot_mat  * R
	bound_R.addVertex(points_R.appendVert( R ))
bound_R.close()



grid.connect(bound_R)
grid.outputPLC( 'out.smesh' )

#pointst = copy.deepcopy(pointsb)

#print points_L.simpleString()
#print points_R.simpleString()

print grid
