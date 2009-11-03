#!/usr/bin/python

from gridhelper import *
import math, copy

area 		= math.pi
tube_length	= 3
num_verts 	= 5
alpha		= math.radians(360 / (num_verts ) )
alpha_half	= alpha / 2.
area_one_tri = math.pi / float(num_verts)
L_x = math.sqrt( area_one_tri*area_one_tri / ( math.sin( alpha_half ) * math.cos( alpha_half ) ) )

origin_L = Vector3(0,0,0)
points_L = PLCPointList( 3 )

bound_L = FanningSimplexList( points_L.appendVert( origin_L ), 3 )
L = Vector3( L_x, 0, 0 )
bound_L.addVertex( points_L.appendVert( L ) )
rot_mat = Matrix4.new_rotatez( alpha )
for i in range( 1, num_verts  ):
	L = rot_mat * L
	bound_L.addVertex(points_L.appendVert( L ))
bound_L.close()

origin_R = Vector3(0,0,tube_length)
points_R = PLCPointList( 3 )
bound_R = FanningSimplexList( points_R.appendVert( origin_R ), 2 )
R = Vector3( L_x, 0, tube_length )
bound_R.addVertex( points_R.appendVert( R ) )
for i in range( 1, num_verts  ):
	R = rot_mat * R
	bound_R.addVertex(points_R.appendVert( R ))
bound_R.close()

grid = FullGrid( bound_L, bound_R, 4 )
grid.connect()
grid.outputPLC( 'out.smesh' )

#pointst = copy.deepcopy(pointsb)

print points_L.simpleString()
print points_R.simpleString()

print 'boundL'
print bound_L
print 'boundR'
print bound_R
print 'grid'
print grid

#out = open( "tube.poly", 'w' )
'''
out.write( '%d 2 0 %d\n'%(nVert,nBids) )
cVert = 0
for v in verts:
		out.write( '%d %f %f 0\n'%(cVert,v[0],v[1]) )
		cVert += 1
out.write( '%d 1\n'%(nBsegs) )
cBseg = 0
for bid,vertIds in bSegs.items():
	for seg in vertIds:
		out.write( '%d %d %d %d\n'%(cBseg,seg[0],seg[1],bid)  )
		cBseg += 1

out.write( '%d\n'%(0))
out.close()
'''